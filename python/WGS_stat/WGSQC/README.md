# 使用说明

pzw

## 说明
程序能统计fastp以及bedtools的下游数据，目的是统计WGS样本的质量。

使用前需手动修改WGSQC_v0.1.py中其他功能脚本的路径，建议使用绝对路径。

### fastp
使用fastp进行样本质控
```bash
fastp -i sample.R1.fq.gz -I sample.R2.fq.gz -o sample.clean.R1.fq.gz -O sample.R2.fq.gz -w 8 -j sample.json -h sample.html
```

程序的fastp功能通过读取sample.json文件统计质量信息
```bash
python WGSQC_v0.1.py fastp -i sample.json -o sample.fastp.txt
```


### bedtools
通过bedtools分析样本的测序深度与覆盖度
**注：比对到hg19.fa的bam文件，理论上染色体表示为“chrx”的样本均可**
```bash
bedtools genomecov -ibam sample.sorted.bam -g hg19.fa -bga > sample.cov.txt
```

程序的bedtools功能可读取sample.cov.txt文件统计覆盖度信息以及平均深度
```bash
python WGSQC_v0.1.py bedtools -i sample.cov.txt -o sample.bedtools.txt
```

也可以单独输出一个染色体的比对信息
```bash
python WGSQC_v0.1.py bedtools -i sample.cov.txt -chr chr1
```

当然，也可以一次指定多个染色体，用“,”分隔。
```bash
python WGSQC_v0.1.py bedtools -i sample.cov.txt -chr chr1,chr3,chr5
```
