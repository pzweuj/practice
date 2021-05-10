# pzw
# 20210510
# hg19
# 用于小panel(组织版)


# configfile:
######## 样本信息，后续改用configfile导入 ########
Sample = "test"
RawdataDir = "test/rawdata"
OutputDir = "test/test"
Threads = 8
################################################


############  数据库  #############
Reference = "/home/bioinfo/ubuntu/database/hg19/ucsc.hg19.fasta"
Mills_indel = "/home/bioinfo/ubuntu/database/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
Indel = "/home/bioinfo/ubuntu/database/hg19/1000G_phase1.indels.hg19.sites.vcf"
Snp = "/home/bioinfo/ubuntu/database/hg19/dbsnp_138.hg19.vcf"
Bed = "/home/bioinfo/ubuntu/database/hg19/cnv/tmerge.bed"
Germline = "/home/bioinfo/ubuntu/database/hg19/af-only-gnomad.raw.sites.hg19.vcf.gz"
Target = "/home/bioinfo/ubuntu/database/hg19/cnv/t.target.bed"
AntiTarget = "/home/bioinfo/ubuntu/database/hg19/cnv/t.antitarget.bed"
NormalCnn = "/home/bioinfo/ubuntu/database/hg19/cnv/reference.t.cnn"
Humandb = "/home/bioinfo/ubuntu/database/humandb"
###################################

#############  软件  ##############
FASTP = "/home/bioinfo/ubuntu/software/fastp/fastp"
BWA = "/home/bioinfo/ubuntu/software/bwa-0.7.17/bwa"
SAMBAMBA = "/home/bioinfo/ubuntu/software/sambamba-0.8.0/sambamba"
GATK = "/home/bioinfo/ubuntu/software/gatk-4.2.0.0/gatk"
CNVKIT = "/home/bioinfo/.local/bin/cnvkit.py"
LUMPYEXPRESS = "/home/bioinfo/ubuntu/software/lumpy-sv/scripts/lumpyexpress"
EXTRACTSPLITREADS = "/home/bioinfo/ubuntu/software/lumpy-sv/scripts/extractSplitReads_BwaMem"
SAMTOOLS = "/home/bioinfo/ubuntu/software/samtools-1.11/samtools"
MANTA = "/home/bioinfo/ubuntu/software/manta-1.6.0/bin/configManta.py"
SNPEFF = "/home/bioinfo/ubuntu/software/snpEff/snpEff.jar"
SNPEFFCONFIG = "/home/bioinfo/ubuntu/software/snpEff/snpEff.config"
COVERT2ANNOVAR = "/home/bioinfo/ubuntu/software/annovar/convert2annovar.pl"
TABLEANNOVAR = "/home/bioinfo/ubuntu/software/annovar/table_annovar.pl"
###################################


# PipeLine
rule all:
    input:
        sv1 = "{OutputDir}/tempFile/lumpy_{Sample}/{Sample}.lumpy.vcf".format(Sample=Sample, OutputDir=OutputDir),
        sv2 = "{OutputDir}/tempFile/manta_{Sample}/{Sample}.manta.vcf".format(Sample=Sample, OutputDir=OutputDir),
        anno = "{OutputDir}/annotation/{Sample}.hg19_multianno.txt".format(Sample=Sample, OutputDir=OutputDir)

# 质控
rule fastp:
    input:
        RawRead1 = "{RawdataDir}/{Sample}_R1.fastq.gz".format(Sample=Sample, RawdataDir=RawdataDir),
        RawRead2 = "{RawdataDir}/{Sample}_R2.fastq.gz".format(Sample=Sample, RawdataDir=RawdataDir)
    output:
        CleanRead1 = "{OutputDir}/cleandata/{Sample}.clean.R1.fastq.gz".format(Sample=Sample, OutputDir=OutputDir),
        CleanRead2 = "{OutputDir}/cleandata/{Sample}.clean.R2.fastq.gz".format(Sample=Sample, OutputDir=OutputDir),
        JsonReport = "{OutputDir}/tempFile/fastp_{Sample}/{Sample}.json".format(Sample=Sample, OutputDir=OutputDir),
        HtmlReport = "{OutputDir}/tempFile/fastp_{Sample}/{Sample}.html".format(Sample=Sample, OutputDir=OutputDir)
    shell:
        """
            mkdir -p {OutputDir}/cleandata
            mkdir -p {OutputDir}/tempFile/fastp_{Sample}
            {FASTP} \\
                -i {input.RawRead1} \\
                -I {input.RawRead2} \\
                -o {output.CleanRead1} \\
                -O {output.CleanRead2} \\
                -j {output.JsonReport} \\
                -h {output.HtmlReport} \\
                -w {Threads}
        """

# 比对
rule bwa:
    input:
        CleanRead1 = rules.fastp.output.CleanRead1,
        CleanRead2 = rules.fastp.output.CleanRead2,
        Reference = "{Reference}".format(Reference=Reference)
    output:
        Bam = "{OutputDir}/tempFile/bwa_{Sample}/{Sample}.bam".format(Sample=Sample, OutputDir=OutputDir)
    shell:
        """
            mkdir -p {OutputDir}
            mkdir -p {OutputDir}/tempFile/bwa_{Sample}
            {BWA} mem -t {Threads} \\
                -R "@RG\\tPL:illumina\\tPU:Test\\tID:{Sample}\\tSM:{Sample}" \\
                {input.Reference} {input.CleanRead1} {input.CleanRead2} \\
                | {SAMBAMBA} view -f bam -t {Threads} -p \\
                -S /dev/stdin > {output.Bam}
        """

# 排序
rule sort:
    input:
        Bam = rules.bwa.output.Bam
    output:
        SortBam = "{OutputDir}/bam/{Sample}.sort.bam".format(Sample=Sample, OutputDir=OutputDir)
    shell:
        """
            mkdir -p {OutputDir}/bam
            {SAMBAMBA} sort {input.Bam} \\
                -t {Threads} -o {output.SortBam} -p
            rm {input.Bam}
        """

# 重复标记
rule markdups:
    input:
        Bam = rules.sort.output.SortBam
    output:
        MarkBam = "{OutputDir}/tempFile/bwa_{Sample}/{Sample}.marked.bam".format(Sample=Sample, OutputDir=OutputDir)
    shell:
        """
            {SAMBAMBA} markdup \\
                {input.Bam} {output.MarkBam} -p --overflow-list-size 600000 \\
                -t {Threads}
        """

# Bam文件质量值校正
rule baseRecal:
    input:
        Bam = rules.markdups.output.MarkBam,
        Snp = "{Snp}".format(Snp=Snp),
        Indel = "{Indel}".format(Indel=Indel),
        Mills_indel = "{Mills_indel}".format(Mills_indel=Mills_indel),
        Reference = "{Reference}".format(Reference=Reference)
    output:
        RecalTable = "{OutputDir}/tempFile/gatk_{Sample}/{Sample}.recal.table".format(Sample=Sample, OutputDir=OutputDir)
    shell:
        """
            mkdir -p {OutputDir}/tempFile/gatk_{Sample}
            {GATK} BaseRecalibrator \\
                --known-sites {input.Snp} \\
                --known-sites {input.Mills_indel} \\
                --known-sites {input.Indel} \\
                -R {input.Reference} \\
                -I {input.Bam} \\
                -O {output.RecalTable}
        """

rule applyBQSR:
    input:
        Bam = rules.markdups.output.MarkBam,
        RecalTable = rules.baseRecal.output.RecalTable,
        Reference = "{Reference}".format(Reference=Reference)
    output:
        RecalBam = "{OutputDir}/bam/{Sample}.Realign.bam".format(Sample=Sample, OutputDir=OutputDir)
    shell:
        """
            {GATK} ApplyBQSR \\
                -R {input.Reference} \\
                --bqsr-recal-file {input.RecalTable} \\
                -I {input.Bam} \\
                -O {output.RecalBam}
            mv {OutputDir}/bam/{Sample}.Realign.bai {OutputDir}/bam/{Sample}.Realign.bam.bai
        """

## 统计质量与比对信息

# SNV/Indel
rule snv_indel:
    input:
        Reference = "{Reference}".format(Reference=Reference),
        Bed = "{Bed}".format(Bed=Bed),
        Bam = rules.applyBQSR.output.RecalBam,
        Germline = "{Germline}".format(Germline=Germline)
    threads: 16
    output:
        Vcf = "{OutputDir}/vcf/{Sample}.vcf".format(Sample=Sample, OutputDir=OutputDir)
    shell:
        """
            mkdir -p {OutputDir}/vcf
            mkdir -p {OutputDir}/tempFile/gatk_{Sample}
            {GATK} Mutect2 \\
                -R {input.Reference} \\
                -I {input.Bam} \\
                -O {OutputDir}/tempFile/gatk_{Sample}/{Sample}.vcf \\
                --germline-resource {input.Germline} \\
                -tumor {Sample} \\
                --native-pair-hmm-threads {Threads} \\
                -L {input.Bed} \\
                -A Coverage -A GenotypeSummaries \\
                --max-reads-per-alignment-start 0
            cp {OutputDir}/tempFile/gatk_{Sample}/{Sample}.vcf {output.Vcf}
        """

# CNV
rule cnv:
    input:
        Bam = rules.applyBQSR.output.RecalBam,
        NormalCnn = "{NormalCnn}".format(NormalCnn=NormalCnn),
        Target = "{Target}".format(Target=Target),
        AntiTarget = "{AntiTarget}".format(AntiTarget=AntiTarget)
    output:
        Cns = "{OutputDir}/tempFile/cnvkit_{Sample}/{Sample}.call.cns".format(Sample=Sample, OutputDir=OutputDir)
    shell:
        """
            {CNVKIT} coverage {input.Bam} \\
                {input.Target} \\
                -o {OutputDir}/tempFile/cnvkit_{Sample}/{Sample}.targetcoverage.cnn \\
                -p {Threads}
            {CNVKIT} coverage {input.Bam} \\
                {input.AntiTarget} \\
                -o {OutputDir}/tempFile/cnvkit_{Sample}/{Sample}.antitargetcoverage.cnn \\
                -p {Threads}
            {CNVKIT} fix \\
                {OutputDir}/tempFile/cnvkit_{Sample}/{Sample}.targetcoverage.cnn \\
                {OutputDir}/tempFile/cnvkit_{Sample}/{Sample}.antitargetcoverage.cnn \\
                {input.NormalCnn} \\
                -o {OutputDir}/tempFile/cnvkit_{Sample}/{Sample}.cnr
            {CNVKIT} segment {OutputDir}/tempFile/cnvkit_{Sample}/{Sample}.cnr \\
                -o {OutputDir}/tempFile/cnvkit_{Sample}/{Sample}.cns -p {Threads}
            {CNVKIT} call {OutputDir}/tempFile/cnvkit_{Sample}/{Sample}.cns \\
                -o {output.Cns}
        """

# SV
rule lumpy:
    input:
        Bam = rules.sort.output.SortBam
    output:
        Vcf = "{OutputDir}/tempFile/lumpy_{Sample}/{Sample}.lumpy.vcf".format(Sample=Sample, OutputDir=OutputDir)
    shell:
        """
            mkdir -p {OutputDir}/tempFile/lumpy_{Sample}
            {SAMTOOLS} view -bh -F 1294 {input.Bam} \\
                | {SAMTOOLS} sort -@ {Threads} - \\
                -o {OutputDir}/tempFile/lumpy_{Sample}/{Sample}.discordants.bam
            {SAMTOOLS} index {OutputDir}/tempFile/lumpy_{Sample}/{Sample}.discordants.bam
            {SAMTOOLS} view -h {input.Bam} \\
                | {EXTRACTSPLITREADS} \\
                -i stdin \\
                | {SAMTOOLS} view -bSh - \\
                | {SAMTOOLS} sort -@ {Threads} - \\
                -o {OutputDir}/tempFile/lumpy_{Sample}/{Sample}.splitters.bam
            {SAMTOOLS} index {OutputDir}/tempFile/lumpy_{Sample}/{Sample}.splitters.bam
            {LUMPYEXPRESS} -B {input.Bam} \\
                -D {OutputDir}/tempFile/lumpy_{Sample}/{Sample}.discordants.bam \\
                -S {OutputDir}/tempFile/lumpy_{Sample}/{Sample}.splitters.bam \\
                -o {output.Vcf}
        """


rule manta:
    input:
        Bam = rules.sort.output.SortBam,
        Reference = "{Reference}".format(Reference=Reference)
    output:
        Vcf = "{OutputDir}/tempFile/manta_{Sample}/{Sample}.manta.vcf".format(Sample=Sample, OutputDir=OutputDir)
    shell:
        """
            mkdir -p {OutputDir}/tempFile/manta_{Sample}
            rm -rf {OutputDir}/tempFile/manta_{Sample}/*
            {MANTA} \\
                --tumorBam {input.Bam} \\
                --referenceFasta {input.Reference} \\
                --exome \\
                --generateEvidenceBam \\
                --runDir {OutputDir}/tempFile/manta_{Sample}
            {OutputDir}/tempFile/manta_{Sample}/runWorkflow.py -j {Threads}
            zcat {OutputDir}/tempFile/manta_{Sample}/results/variants/tumorSV.vcf.gz \\
                > {OutputDir}/tempFile/manta_{Sample}/{Sample}.manta.vcf
        """

# 注释
rule annotation:
    input:
        Vcf = rules.snv_indel.output.Vcf,
        Humandb = "{Humandb}".format(Humandb=Humandb)
    output:
        SnpEffResults = "{OutputDir}/annotation/{Sample}.vcf".format(Sample=Sample, OutputDir=OutputDir),
        AnnovarResults = "{OutputDir}/annotation/{Sample}.hg19_multianno.txt".format(Sample=Sample, OutputDir=OutputDir)
    shell:
        """
            mkdir -p {OutputDir}/tempFile/anno_{Sample}
            java -jar {SNPEFF} -c {SNPEFFCONFIG} \\
                -s {OutputDir}/tempFile/anno_{Sample}/{Sample}.summary.html \\
                hg19 {input.Vcf} > {output.SnpEffResults}
            {COVERT2ANNOVAR} -format vcf4 \\
                {output.SnpEffResults} --includeinfo > {OutputDir}/tempFile/anno_{Sample}/{Sample}.avinput
            {TABLEANNOVAR} {OutputDir}/tempFile/anno_{Sample}/{Sample}.avinput \\
                {input.Humandb} -buildver hg19 -out {OutputDir}/tempFile/anno_{Sample}/{Sample} \\
                -remove \\
                -protocol refGene,avsnp150,gnomad211_genome,clinvar_20210308,JaxCkb,Civic,OncoKB,dbnsfp41a,cosmic92_coding,intervar_20180118 \\
                -operation g,f,f,f,f,f,f,f,f,f \\
                -nastring - -thread {Threads} -otherinfo
            cp {OutputDir}/tempFile/anno_{Sample}/{Sample}.hg19_multianno.txt {output.AnnovarResults}
        """