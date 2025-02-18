# coding=utf-8
# pzw
# 20250218

import streamlit as st
import os
import re
import logging
from logging.handlers import RotatingFileHandler
from pdf2markdown_core import pdf2markdown

# 配置日志记录
def setup_logger():
    formatter = logging.Formatter(
        fmt='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    file_handler = RotatingFileHandler(
        'app.log',
        maxBytes=10 * 1024 * 1024,
        backupCount=2,
        encoding='utf-8'
    )
    file_handler.setFormatter(formatter)
    logger = logging.getLogger('pdf_converter')
    logger.addHandler(file_handler)
    logger.setLevel(logging.INFO)
    return logger

logger = setup_logger()

# 清理文件名中的无效字符
def sanitize_filename(filename):
    invalid_chars = r'["\\/|?:*<>™{}[\]()「」【】『』·]'
    cleaned_name = re.sub(r'[{}]'.format(invalid_chars), '', filename)
    return cleaned_name if cleaned_name else 'untitled'

def main():
    st.set_page_config(
        page_title="PDF转Markdown工具",
        page_icon="📄",
        layout="centered",
        initial_sidebar_state="collapsed"
    )

    st.markdown("""
    <style>
        .stProgress > div > div > div > div { background-color: #4CAF50; }
        .st-bb { background-color: #f0f2f6; }
        .st-at { background-color: #ffffff; }
        .css-18e3th9 { padding: 2rem 1rem; }
    </style>
    """, unsafe_allow_html=True)

    st.title("📄 PDF转Markdown工具")
    st.caption("使用阿里百炼API将PDF文档转换为结构化的Markdown格式")

    # 创建固定工作目录
    working_dir = os.path.join(os.getcwd(), ".temp")
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
        logger.info(f"创建工作目录: {working_dir}")

    with st.form("conversion_form"):
        uploaded_file = st.file_uploader(
            "上传PDF文件",
            type=["pdf"],
            help="请选择需要转换的PDF文档"
        )
        col1, col2 = st.columns(2)
        with col1:
            title = st.text_input(
                "文档标题",
                value="默认标题",
                help="设置输出Markdown的一级标题"
            )
            api_key = st.text_input(
                "API密钥",
                type="password",
                help="阿里百炼的API访问密钥"
            )
        with col2:
            time_gap = st.number_input(
                "请求间隔（秒）",
                min_value=1,
                value=2,
                help="API请求之间的最小间隔时间"
            )
            keep_temp = st.checkbox(
                "保留临时文件",
                help="保留转换过程中生成的临时文件",
                value=False
            )
        if submitted := st.form_submit_button("开始转换", use_container_width=True):
            if not uploaded_file or not api_key:
                st.error("请确保已上传PDF文件并填写API密钥。")
                return

            try:
                with st.spinner("初始化转换环境..."):
                    # 保存上传的文件到指定目录
                    original_filename = uploaded_file.name
                    sanitized_filename = sanitize_filename(original_filename)
                    
                    st.write(f"原始文件名: {original_filename}")
                    st.write(f"清理后的文件名: {sanitized_filename}")
                    
                    input_path = os.path.join(working_dir, sanitized_filename)
                    logger.info(f"文件保存路径: {input_path}")
                    
                    with open(input_path, "wb") as f:
                        f.write(uploaded_file.getbuffer())
                    
                    if not os.path.exists(input_path):
                        logger.error(f"文件 {input_path} 未找到")
                        st.error(f"保存文件失败，路径: {input_path}")
                        return

                    # 设置输出路径
                    output_path = os.path.join(working_dir, "output.md")
                    logger.info(f"输出路径: {output_path}")

                    # 使用try-except包围pdf2markdown的调用
                    try:
                        with st.spinner("进行PDF到Markdown转换..."):
                            # 调用pdf2markdown函数
                            pdf2markdown(
                                input_pdf=input_path,
                                output_md=output_path,
                                title=title,
                                api_key=api_key,
                                time_gap=time_gap,
                                tmp=keep_temp
                            )
                            logger.info("PDF转Markdown完成")

                        # 验证输出文件
                        if os.path.exists(output_path):
                            logger.info(f"输出文件已创建: {output_path}")
                            with open(output_path, "r", encoding="utf-8") as f:
                                md_content = f.read()
                            st.success("转换完成。请下载结果文件。")
                            st.download_button(
                                label="下载Markdown",
                                data=md_content,
                                file_name=f"{sanitized_filename}_converted.md",
                                mime="text/markdown",
                                use_container_width=True
                            )
                        else:
                            logger.error(f"输出文件 {output_path} 未创建")
                            st.error(f"转换失败: 输出文件 {output_path} 未找到")
                    except Exception as e:
                        logger.error(f"调用pdf2markdown失败: {str(e)}", exc_info=True)
                        st.error(f"转换失败: {str(e)}")

            except Exception as e:
                logger.error(f"整体流程失败: {str(e)}", exc_info=True)
                st.error(f"转换失败: {str(e)}")

if __name__ == "__main__":
    main()

# end
