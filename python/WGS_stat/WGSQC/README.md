# 统计WGS质控与比对信息



### 质控信息
质控信息使用原始输入为fastp
```bash
fastp -i sample.R1.fq.gz -I sample.R2.fq.gz -o sample.clean.R1.fq.gz -O sample.clean.R2.fq.gz -w 8 -j sample.json -h sample.html
```
将sample.json作为原始输入。

### 比对信息
可输出比对到各个染色体的覆盖度与平均深度，原始输入为bedtools结果
```bash
bedtools genomecov -ibam sample.bam -g hg19.fa -bga > sample.cov.txt
```


### 使用说明

```bash
python WGSQC_v0.1.py -h
```