# coding=utf-8
# pzw
# 20230620
# 文本加密与解密
# 加密：把文本文件加密后，生成一个新的文本文件
# 将byte使用utf-8编码
# 基于AES-256-GCM
# 解密：1 把加密的文本文件解密，生成解密后的文本文件；2 把加密的文本文件读入，并解密为一个list

import os
import sys
import base64
import hashlib
from Crypto.Cipher import AES
from Crypto.Util.Padding import pad
from Crypto import Random

################################# AES 256 算法 ######################################
# 字符串加密
def encrypt_AES(raw, password):
    private_key = hashlib.sha256(password.encode("utf-8")).digest()
    raw = pad(raw.encode("utf-8"), 16, style="pkcs7")
    iv = Random.new().read(AES.block_size)
    cipher = AES.new(private_key, AES.MODE_GCM, iv)
    encrypt_data = base64.b64encode(iv + cipher.encrypt(raw))
    return encrypt_data.decode("utf-8")

# 字符串解密
def decrypt_AES(enc, password):
    enc = enc.encode("utf-8")
    unpad = lambda x: x[: -ord(x[len(x) - 1:])]
    private_key = hashlib.sha256(password.encode("utf-8")).digest()
    enc = base64.b64decode(enc)
    iv = enc[:16]
    cipher = AES.new(private_key, AES.MODE_GCM, iv)
    unppad = unpad(cipher.decrypt(enc[16:]))
    return unppad.decode("utf-8")

# 文本文件加密
def file_encrypt_AES(inputFile, outputFile, password):
    input = open(inputFile, "r", encoding="utf-8")
    output = open(outputFile, "w", encoding="utf-8")
    
    for line in input:
        encryptString = encrypt_AES(line.rstrip(), password)
        output.write(encryptString + "\n")
    
    input.close()
    output.close()

# 文本文件解密
def file_decrypt_AES(inputFile, outputFile, password):
    input = open(inputFile, "r", encoding="utf-8")
    output = open(outputFile, "w", encoding="utf-8")

    for line in input:
        outputString = decrypt_AES(line.rstrip(), password)
        output.write(outputString + "\n")

    output.close()
    input.close()
    
# 读取加密的文本文件
def file_decrypt_read_AES(inputFile, password):
    input = open(inputFile, "r", encoding="utf-8")
    lineList = []
    for line in input:
        outputString = decrypt_AES(line.strip(), password)
        lineList.append(outputString)
    input.close()
    return lineList

# Main
def main(inputFileString, password, mode):
    # 判断文件是否存在，不存在则认为是字符串
    input = ""
    if os.path.exists(inputFileString):
        input = "File"
    else:
        input = "String"

    if mode == "E":
        if input == "String":
            outputString = encrypt_AES(inputFileString, password)
            print(outputString)
        elif input == "File":
            file_encrypt_AES(inputFileString, inputFileString + ".encrypt", password)
        else:
            print("[ERROR] Unknown Input Type!")
    elif mode == "D":
        if input == "String":
            try:
                outputString = decrypt_AES(inputFileString, password)
                print(outputString)
            except:
                print("[ERROR] Decrypt Failed!")
        elif input == "File":
            try:
                file_decrypt_AES(inputFileString, inputFileString + ".decrypt", password)
            except:
                print("[ERROR] Decrypt Failed!")
        else:
            print("[ERROR] Unknown Input Type!")
    else:
        print("[ERROR] Unknown Mode: E/D")

##########
try:
    main(sys.argv[1], sys.argv[2], sys.argv[3])
except:
    print("Usage: python3 EnDeC.py <input> <password> <mode>")
