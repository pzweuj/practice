# coding=utf-8
# pzw
# 20250218
# 使用阿里云百炼的API
# 使用模型 qwen2.5-vl-7b-instruct
# 实测7B模型效果可接受，如需提升效果，可使用72B模型

import os
import sys
import fitz
from PIL import Image
import base64
from openai import OpenAI, APIError
import time
from datetime import datetime
import re
import argparse
import shutil

# PDF 2 Image
def pdf_to_images(pdf_path, output_folder):
    """
    将 PDF 文件的每一页保存为 JPEG 图像，并存储到指定文件夹中。
    
    :param pdf_path: PDF 文件路径
    :param output_folder: 输出图像的文件夹路径
    """
    # 确保输出文件夹存在
    os.makedirs(output_folder, exist_ok=True)
    
    # 打开 PDF 文件
    pdf_document = fitz.open(pdf_path)
    print(f"[Process] Total pages: {pdf_document.page_count}")
    
    # 遍历每一页并将其转换为图像
    for page_num in range(pdf_document.page_count):
        # 获取当前页
        page = pdf_document.load_page(page_num)
        
        # 将页面渲染为图像 (Pixmap)
        pix = page.get_pixmap()
        
        # 将 Pixmap 转换为 PIL 图像
        img = Image.frombytes("RGB", [pix.width, pix.height], pix.samples)
        
        # 构造输出文件名
        image_filename = os.path.join(output_folder, f"page_{page_num + 1}.jpeg")
        
        # 保存图像为 JPEG 格式
        img.save(image_filename, "JPEG")
        print("[Process] Page: " +  f"page_{page_num + 1}")
    print("[Process] PDF 2 Images Done")

# 阿里百炼API调用
def qwen_vl_api(input_img, api_key):
    model = "qwen2.5-vl-7b-instruct"
    api_base_url = "https://dashscope.aliyuncs.com/compatible-mode/v1"
    
    # 获取图片路径并读取为 Base64 编码
    input_img_path = os.path.abspath(input_img)
    with open(input_img_path, "rb") as image_file:
        base64_image = base64.b64encode(image_file.read()).decode("utf-8")
    
    client = OpenAI(
        api_key=api_key,
        base_url=api_base_url,
    )
    
    prompt = """
    请仔细分析图片内容，将图片中的文字和表格信息转换为规范的markdown格式。具体要求如下：
    1. 识别所有文字内容，按原格式保留段落、标题、列表等结构
    2. 识别表格内容，使用markdown表格语法进行转换
    3. 识别数学公式，使用markdown公式语法进行转换
    4. 保留原始文本的层次结构，使用#、##等标记表示标题级别
    5. 对于列表内容，使用-或*标记
    6. 确保转换后的markdown格式规范，可直接用于文档编写
    7. 如果图片中包含代码块，请使用```标记进行包裹
    8. 非代码块内容，不要用```标记进行包裹
    9. 如果是一张不包含信息的图片，则使用<img_start> <img_end>标签进行标记，并在标签内描述图片的内容
    10. 保持原始内容的准确性，不要添加或删除任何信息
    """
    
    max_retries = 5  # 最大重试次数
    retry_count = 0  # 当前重试次数
    
    while retry_count < max_retries:
        try:
            completion = client.chat.completions.create(
                model=model,
                messages=[{"role": "user", "content": [
                    {"type": "text", "text": prompt},
                    {"type": "image_url",
                     "image_url": {"url": f"data:image/jpeg;base64,{base64_image}"}}
                ]}]
            )
            
            # 解析返回结果
            markdown_content = completion.choices[0].message.content
            
            # 将 Markdown 内容封装为一个字典对象
            markdown_object = {
                "title": "Extracted Content from Image",
                "content": markdown_content,
                "metadata": {
                    "model": model,
                    "input_image_path": "file://" + input_img_path,
                    "api_base_url": api_base_url
                }
            }
            return markdown_object
        
        except APIError as e:
            retry_count += 1
            print(f"请求失败，正在进行第 {retry_count} 次重试。错误信息: {e}")
            if retry_count >= max_retries:
                print("[Error] ", input_img_path, "已达到最大重试次数，请求仍然失败。")
                raise  # 抛出异常以终止程序

# 保存为一个markdown
def save_to_md(content, output_file):
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(content)
        f.write("\n\n")

# 批量转换
def all_images_to_md(input_dir, output_dir, api_key, time_gap=2):
    os.makedirs(output_dir, exist_ok=True)
    n = 1
    while n > 0:
        img = os.path.join(input_dir, f"page_{n}.jpeg")
        output = os.path.join(output_dir, f"page_{n}.md")
        if os.path.exists(img):
            result = qwen_vl_api(img, api_key)
            md_content = result["content"]
            save_to_md(md_content, output)
            print(f"[Process] Page {n} VL Done")
        else:
            break
        n += 1
        time.sleep(time_gap)

# 合并一个文件夹下的所有markdown文件
def merge_all_markdown(input_dir, output_md, title="默认标题"):
    with open(output_md, "w", encoding="utf-8") as f:
        f.write("# " + title.lstrip("#") + "\n")
        n = 1
        while n > 0:
            page = os.path.join(input_dir, f"page_{n}.md")
            if os.path.exists(page):
                with open(page, "r", encoding="utf-8") as p:
                    for line in p:
                        line = re.sub(r'^(#+) ', r'#\1 ', line)
                        f.write(line)
            else:
                break
            n += 1

# 主函数
def pdf2markdown(input_pdf, output_md, title, api_key, time_gap, tmp):
    input_pdf = os.path.abspath(input_pdf)
    output_md = os.path.abspath(output_md)
    output_dir = os.path.dirname(output_md)
    date = datetime.now().strftime("%d-%m-%Y-%H-%M-%S")
    image_tmp_path = os.path.join(output_dir, date + "_img_tmp")
    md_tmp_path = os.path.join(output_dir, date + "_md_tmp")
    pdf_to_images(input_pdf, image_tmp_path)
    all_images_to_md(image_tmp_path, md_tmp_path, api_key, time_gap)
    merge_all_markdown(md_tmp_path, output_md, title)

    # 删除临时文件夹
    if not tmp:
        shutil.rmtree(image_tmp_path)
        shutil.rmtree(md_tmp_path)

    print("[Done]", output_md)

def main():
    parser = argparse.ArgumentParser(
        description="使用示例：python3 pdf2markdown_core.py -i <input pdf> -o <output markdown> -t <title> -k <api key>",
        prog="pdf2markdown_core.py",
        usage="python3 pdf2markdown_core.py [-h] -i <input pdf> -o <output markdown> -t <title> -k <api key>",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("-v", "--version", action="version", version="Version 0.1 20250214")
    parser.add_argument("-i", "--input", type=str, help="输入pdf")
    parser.add_argument("-o", "--output", type=str, help="输出markdown")
    parser.add_argument("-t", "--title", type=str, help="输出markdown文件的一级标题", default="默认标题")
    parser.add_argument("-k", "--key", type=str, help="阿里百炼API KEY")
    parser.add_argument("-g", "--gap", type=int, help="每个对话的时间间隔秒数，默认是2秒", default=2)
    parser.add_argument("-r", "--tmp", type=bool, help="是否保留临时文件夹，默认是False [True/False]", default=False)
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()
    pdf2markdown(input_pdf=args.input, output_md=args.output, title=args.title, api_key=args.key, time_gap=args.gap, tmp=args.tmp)

if __name__ == "__main__":
    main()

# end
