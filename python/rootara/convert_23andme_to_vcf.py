# coding=utf-8
# pzw
# 20250114

import sys
import pandas as pd
import pysam
import argparse

def convert_23andme_to_vcf(input_file, output_file, reference_fasta):
    # 读取23andMe数据文件
    df = pd.read_csv(input_file, sep='\t', comment='#', header=None, names=['rsid', 'chromosome', 'position', 'genotype'], low_memory=False)

    # 过滤掉无效的行（例如，染色体为0的行）
    df = df[df['chromosome'] != '0']

    # 将染色体名称转换为VCF格式（例如，'23' -> 'X', '24' -> 'Y', '25' -> 'MT'）
    df['chromosome'] = df['chromosome'].replace({'23': 'X', '24': 'Y', '25': 'MT'})

    # 打开参考基因组FASTA文件
    fasta = pysam.FastaFile(reference_fasta)

    # 创建VCF文件头
    vcf_header = f"""##fileformat=VCFv4.2
##source=23andMe
##reference={reference_fasta}
##INFO=<ID=RS,Number=1,Type=String,Description="dbSNP ID">
##INFO=<ID=RAW,Number=1,Type=String,Description="Raw Genotype">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
"""

    # 打开输出文件并写入VCF头
    with open(output_file, 'w') as vcf_file:
        vcf_file.write(vcf_header)

        # 遍历每一行数据并写入VCF格式
        for _, row in df.iterrows():
            chrom = row['chromosome']
            pos = int(row['position'])
            rsid = row['rsid']
            genotype = row['genotype']

            if genotype == "--":
                continue

            # 从参考基因组中获取参考等位基因
            try:
                # 确保chromosome是字符串类型
                chrom_str = str(chrom)
                # 获取参考等位基因，pysam使用0-based坐标
                ref_allele = fasta.fetch(chrom_str, pos - 1, pos).upper()
            except KeyError:
                print(f"Warning: Chromosome {chrom} not found in reference genome. Skipping.")
                continue

            # 确定ALT等位基因
            alt_alleles = set(genotype) - {ref_allele}
            alt_allele = ','.join(alt_alleles) if alt_alleles else '.'

            # 确定基因型GT
            if chrom in ['MT', 'Y']:  # 处理线粒体和Y染色体
                if alt_allele == '.':
                    gt = '0'
                else:
                    gt = '1'
            else:  # 处理常染色体和X染色体
                if alt_allele == '.':  # 野生型
                    gt = '0/0'
                elif len(genotype) == 1:  # 纯合突变
                    gt = '1/1'
                else:  # 杂合突变
                    if len(set(genotype)) == 1:  # 纯合
                        gt = '1/1'
                    else:  # 杂合
                        gt = '0/1'
            
            # 修复alt
            if gt in ['0/0', '0']:
                alt_allele = ref_allele

            # 写入VCF行
            vcf_file.write(f"{chrom}\t{pos}\t{rsid}\t{ref_allele}\t{alt_allele}\t.\t.\tRS={rsid};RAW={genotype}\tGT\t{gt}\n")

    # 关闭参考基因组文件
    fasta.close()

if __name__ == "__main__":
    # 创建参数解析器
    parser = argparse.ArgumentParser(description='Convert 23andme To Vcf')
    
    # 添加参数
    parser.add_argument('-i', '--input', required=True, help='23andMe Data')
    parser.add_argument('-o', '--output', required=True, help='output.vcf')
    parser.add_argument('-r', '--reference', required=True, help='reference.fasta')
    
    # 如果没有传递任何参数，打印帮助文档
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    
    # 解析参数
    args = parser.parse_args()
    
    # 调用转换函数
    convert_23andme_to_vcf(args.input, args.output, args.reference)
    print(f"VCF Created: {args.output}")

# end
