# coding=utf-8
# pzw
# 20230103
# 扫描目录是否下机完成，并把数据推送到分析服务器中

import os
import time
from datetime import datetime

# 获得文件夹大小
def getDirSize(dir):
    size = 0
    for root, dirs, files in os.walk(dir):
        size += sum([os.path.getsize(os.path.join(root, name)) for name in files])
    return size

# 已下机的批次
def main():
    tNGSFinished = []
    with open("NGS.txt", "r", encoding="utf-8") as t:
        for line in t:
            tNGSFinished.append(line.replace("\n", ""))

    now = datetime.now()

    path = "/Uploads"
    fileList = os.listdir(path)
    for f in fileList:
        try:
            chipName = f.split("_")[1]
        except:
            chipName = "-"
        # 确认不在已经下机的批次中
        if not chipName in tNGSFinished:
            # 判断是否拆分
            # basecall/A1201021797_OutputFq/mergeDies
            if os.path.exists("sending.txt"):
                continue
            elif os.path.exists(path + "/" + f + "/basecall/" + chipName + "_OutputFq/mergeDies"):
                outputPath = path + "/" + f + "/basecall/" + chipName + "_OutputFq/mergeDies"
            elif os.path.exists(path + "/" + f + "/basecall/OutputFq/mergeDies"):
                outputPath = path + "/" + f + "/basecall/OutputFq/mergeDies"
            else:
                continue
            
            # 确认是文件夹
            if os.path.isdir(outputPath):
                # 确认文件夹的大小
                s = getDirSize(outputPath)
                # 记录下当前的大小
                if "size.txt" in os.listdir("/autoSCP"):
                    sizeFile = open("size.txt", "r", encoding="utf-8")
                    chipNameCheck = "-"
                    checkSize = "10"
                    for line in sizeFile:
                        lines = line.replace("\n", "").split("\t")
                        chipNameCheck = lines[0]
                        checkSize = lines[1]
                    sizeFile.close()
                    # 就是已经跑过，正在下机
                    if chipName == chipNameCheck:
                        # 大小不变了，认为已经下机完
                        if str(s) == checkSize:
                            # 创建一个sending文件，避免同时传输
                            sending = open("sending.txt", "w")
                            sending.close()

                            cmd = """
                                ssh user@192.168.1.1 "mkdir -p /output/{chipName}/L01"
                            """.format(chipName=chipName)
                            os.system(cmd)
                            fqList = os.listdir(outputPath)
                            for fq in fqList:
                                if "_Barcode_" in fq:
                                    rename = fq.replace("_Barcode_", "_L01_")
                                    cmd = """
                                        scp {outputPath}/{fq} user@192.168.1.1:/output/{chipName}/L01/{rename}
                                    """.format(outputPath=outputPath, fq=fq, chipName=chipName, rename=rename)
                                    os.system(cmd)
                                else:
                                    cmd = """
                                        scp {outputPath}/{fq} user@192.168.1.1:/output/{chipName}/L01/
                                    """.format(outputPath=outputPath, fq=fq, chipName=chipName)
                                    os.system(cmd)
                            cmd = """
                                ssh user@192.168.1.1 "touch /output/Info/Upload/{chipName}_L01_DualBarcode_Success.txt"
                            """.format(chipName=chipName)
                            os.system(cmd)
                            print("【完成】批次" + chipName + "推送完成", now)
                            os.remove("sending.txt")
                            with open("NGS.txt", "a", encoding="utf-8") as t:
                                t.write(chipName + "\n")
                        # 大小还在变，就更新大小
                        else:
                            sizeFile = open("size.txt", "w", encoding="utf-8")
                            sizeFile.write(chipName + "\t" + str(s) + "\n")
                            sizeFile.close()
                    # 新的下机批次
                    else:
                        print("【下机】检测到新批次正在下机" + chipName, now)
                        sizeFile = open("size.txt", "w", encoding="utf-8")
                        sizeFile.write(chipName + "\t" + str(s) + "\n")
                        sizeFile.close()
                else:
                    print("【下机】检测到新批次正在下机" + chipName, now)
                    sizeFile = open("size.txt", "w", encoding="utf-8")
                    sizeFile.write(chipName + "\t" + str(s) + "\n")
                    sizeFile.close()     

# 测序仪是离线的，装不了MTA导致crontab无法使用，这里直接使用python来定时
while True:
    main()
    # 5分钟执行一次
    time.sleep(5*60)

# end


