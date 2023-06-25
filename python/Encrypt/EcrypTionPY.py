#!/huayin/software/INDEX/python3
# coding=utf-8
# pzw
# 20230530
# 基于python3
# 使用cython对python脚本进行加密
# 在原脚本中必须将模块写为函数
# 然后在运行脚本中import进需求的模块（注意运行脚本不进行加密，最后使用Py-Fuscate进行混淆）
# 运行脚本应当尽量简洁
# 总体架构如下
#######################################################
# - function\
# - - module1.py
# - - module2.py
# - - module3.py
# - run.py

##
# 在run.py中
# import function.module1.definition1
# import function.module2.definition2
# import function.module3.definition3

## 
# 使用方式
# python EcrypTionPY.py -i <input_script.py> -o <output_directory>
# 最终会在<output_directory>生成<input_script>文件，可使用
# python3 <input_script>
# 建议在<input_script.py>头中指定
# #!/usr/bin/env python3
# 这样在chmod +x <input_script>后
# 可直接使用<input_script>来调用
#######################################################

import os
import sys
import shutil
import argparse

# 获得脚本所在路径
def getCurrentPath(filename):
    now = os.path.abspath(os.path.dirname(filename))
    return now

# 获得程序所在路径
def getAbsPath():
    now = os.path.abspath(os.path.dirname(sys.argv[0]))
    return now

# cython编译脚本生成
def cythonCompile(script, outputDir):
    # 当前系统环境的后缀
    endName = ".cpython-37m-x86_64-linux-gnu.so"

    baseName = os.path.basename(script).replace(".py", "")
    c_script = script.replace(".py", ".c")
    setup_script = open(outputDir + "/" + baseName + ".setup.py", "w", encoding="utf-8")
    setup_script.write("from distutils.core import setup\n")
    setup_script.write("from Cython.Build import cythonize\n")
    setup_script.write("module = cythonize(['{s}'], build_dir='{d}', language_level='3')\n".format(s=script.replace("\\", "/"), d=outputDir.replace("\\", "/")))
    setup_script.write("setup(ext_modules = module)\n")
    setup_script.close()
    try:
        os.remove(outputDir + "/" + baseName + endName)
    except:
        pass
    os.system("python3 " + outputDir + "/" + baseName + ".setup.py build_ext --inplace")
    
    os.remove(outputDir + "/" + baseName + ".setup.py")
    os.remove(c_script)
    # os.remove("build")
    
    # cython会在当前cmd目录下生成pyd，要将对应的pyd移动回去
    shutil.move(baseName + endName, outputDir + "/")

# 要加密run.py
def findImportListAndFuscate(inputScript, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    ready_encrypt = open(inputScript, "r", encoding="utf-8")
    importList = []
    for line in ready_encrypt:
        if line.strip().startswith("#"):
            continue
        elif line.strip().startswith("import"):
            importList.append(line.strip().replace("import ", ""))
        elif line.strip().startswith("from"):
            importList.append(line.replace("from ", "").split(" import ")[0])
        else:
            continue

    ready_encrypt.close()
    currentPath = getCurrentPath(inputScript)
    g = os.walk(currentPath)

    for path, dir_list, file_list in g:
        for file_name in file_list:
            if file_name.endswith(".py"):
                # print(file_name)
                file_path = os.path.join(path, file_name)
                # 这是引用的脚本
                if file_path.replace(currentPath + "/", "").replace(".py", "").replace("/", ".") in importList:
                    script_path = getCurrentPath(file_path)
                    if script_path == currentPath:
                        cythonCompile(file_path, output_dir)
                    else:
                        output_path = output_dir + script_path.replace(currentPath, "")
                        if not os.path.exists(output_path):
                            os.makedirs(output_path, exist_ok=True)
                        cythonCompile(file_path, output_path)
                    findImportListAndFuscate(file_path, script_path)
    
# 加密流程
def main(inputScript, outputDir):
    # 必须是python文件
    if not inputScript.endswith(".py"):
        print("【请输入python文件】")
        exit()

    # 首先对原始输入的script中import的库进行加密
    findImportListAndFuscate(inputScript, outputDir)

    # 然后再对这个inputScript进行混淆
    filename = os.path.basename(inputScript).replace(".py", "")
    output_path = outputDir + "/" + filename + ".py"
    fuscatePath = getAbsPath() + "/Py-Fuscate-main/py_fuscate.py"
    cmd = """
        python3 {fuscatePath} -i {inputScript} -o {output_path} -c 50
    """.format(inputScript=inputScript, output_path=output_path, fuscatePath=fuscatePath)
    os.system(cmd)

# Main
if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description="EnCrypTion Python Script @PZW",
		prog="EcrypTionPY.py",
		usage="python3 EcrypTionPY.py [-h] -i <python3 script> -o <output directory>",
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser.add_argument("-v", "--version", action="version",
		version="Version 0.1 20230531")
	parser.add_argument("-i", "--script", type=str,
		help="Input the python script with entension '.py'")
	parser.add_argument("-o", "--output", type=str,
		help="Output directory")

	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(inputScript=args.script, outputDir=args.output)

# end
