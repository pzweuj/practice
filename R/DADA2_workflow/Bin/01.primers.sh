for sample in $(cat ../Rawdata/samples.txt); \
do bbduk.sh \
	in=../Rawdata/"$sample"_sub_R1.fq \
	in2=../Rawdata/"$sample"_sub_R2.fq \
	out=../Filtdata/"$sample"_sub_R1_trimmed.fq.gz \
	out2=../Filtdata/"$sample"_sub_R2_trimmed.fq.gz \
	literal=GTGCCAGCMGCCGCGGTAA,GGACTACHVGGGTWTCTAAT \
	k=10 \
	ordered=t \
	mink=2 \
	ktrim=l \
	rcomp=f \
	minlength=220 \
	maxlength=280 \
	tbo=t \
	tpe=t; \
done 2> ../Filtdata/bbduk_primer_trim.txt