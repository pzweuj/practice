##
## usage : sh fastq2fasta.sh in.fastq out.fasta

cat $1 | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' > $2