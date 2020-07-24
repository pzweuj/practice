# coding=utf-8
# pzw
# 20200724

from PyPDF2 import PdfFileWriter
from PyPDF2 import PdfFileReader
from PyPDF2 import PdfFileMerger

# 删除pdf中的某些页
def delete_pdf_page(input, output, page_list):
	inputPdf = open(input, "rb")
	inputFile = PdfFileReader(inputPdf)
	outputPdf = open(output, "wb")
	outputFile = PdfFileWriter()
	pages = inputFile.getNumPages()
	pl = page_list.split(",")

	for i in range(pages):
		if str(i + 1) not in pl:
			outputFile.addPage(inputFile.getPage(i))

	outputFile.write(outputPdf)
	inputPdf.close()
	outputPdf.close()

# 用法：
# delete_pdf_page("ABO.pdf", "ABO.new.pdf", "2,3")


# 按顺序合并pdf
def merge_pdf(pdfList, output):
	merger = PdfFileMerger()

	pl = pdfList.split(",")
	for p in pl:
		print(p)
		x = open(p, "rb")
		merger.append(fileobj=x) # 加入到pdf中
		# merger.merge(position=2, fileobj=x) # 从2号位置开始插入
		# merger.appedn(fileobj=x, pages=(0, 2)) # 加入前3页

	outputPdf = open(output, "wb")
	merger.write(outputPdf)
	outputPdf.close()

# 用法
# merge_pdf("ABO.pdf,ABO.new.pdf,ABO.merge.pdf", "ABO.merge2.pdf")

# 顺时针选择90度
def rotate90(input, index):
	inputPdf = open(input, "rb")
	inputFile = PdfFileReader(inputPdf)
	inputFile.getPages(index).rotateClockwise(90)