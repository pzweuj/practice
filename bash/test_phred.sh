##
## usage: sh test_phred in.fastq

less $1 | head -n 1000 | awk '{if(NR%4==0) printf("%s",$0);}' \
	| od -A n -t u1 -v \
	| awk 'BEGIN{min=100;max=0;} \
	{for(i=1;i<=NF;i++) {if($i>max) max=$i; if($i<min) min=$i;}}END \
		{if(max<=126 && min<59) print "Phred33"; \
		else if(max>73 && min>=64) print "Phred64"; \
		else if(min>=59 && min<64 && max>73) print "Solexa64"; \
		else print "Unknown score encoding"; \
			print "( " min ", " max, ")";}'
