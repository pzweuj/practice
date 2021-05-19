version 1.0
# pzw
# 20210519
# Annotation


# snpEff
# https://pcingola.github.io/SnpEff/
task Snpeff {
    input {
        String sample
        File vcf
    }

    File SNPEFF = "/home/bioinfo/ubuntu/software/snpEff/snpEff.jar"
    File config = "/home/bioinfo/ubuntu/software/snpEff/snpEff.config"

    command <<<
        java -jar ~{SNPEFF} \
            -c ~{config} \
            -s ~{sample}.summary.html \
            hg19 ~{vcf} > ~{sample}.snpeff.vcf
    >>>

    output {
        File annoVcf = "~{sample}.snpeff.vcf"
        File summary = "~{sample}.summary.html"
    }

}

# annovar
# https://annovar.openbioinformatics.org/en/latest/
task Annovar {
    input {
        String sample
        File vcf
        Int threads
    }

    String humandb = "/home/bioinfo/ubuntu/database/humandb"

    command <<<
        convert2annovar.pl -format vcf4 \
            ~{vcf} \
            --includeinfo > ~{sample}.avinput
        table_annovar.pl ~{sample}.avinput \
            ~{humandb} -buildver hg19 \
            -out ~{sample} -remove \
            -protocol refGene,avsnp150,gnomad211_genome,clinvar_20210308,JaxCkb,Civic,OncoKB,dbnsfp41a,cosmic92_coding,intervar_20180118 \
            -operation g,f,f,f,f,f,f,f,f,f \
            -nastring - -thread ~{threads} -otherinfo
    >>>

    output {
        File annovarResults = "~{sample}.hg19_multianno.txt"
    }
}
