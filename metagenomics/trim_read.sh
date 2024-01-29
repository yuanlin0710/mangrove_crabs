DIR=/home/lynn/data/Tom_fyp_gut_meta
adaptor=${DIR}/adapter3.fa
for sample in $DIR/*_R1.fq.gz
do
name=${sample%_*1*}
var=$(echo $name |cut -d "/" -f6)
echo $sample $var ${name}_R2.fq.gz
fastp -i $sample -I ${name}_R2.fq.gz -o $var.trim.R1.fq.gz -O $var.trim.R2.fq.gz --adapter_fasta $adaptor -w 64 -h ${outname}_fastp.html --correction
done
