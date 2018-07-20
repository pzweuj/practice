##########
inputdata=NPC10F-N_1.fastq
##########

awk '{if(NR%4==2) print length($1)}' $inputdata | sort -n | uniq -c > reads_length.txt
