awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }' ./${target_name}.all.maker.noseq.gff | \
awk -v OFS="\t" '{ print $1, $2, $3 }' | \
bedtools getfasta -fi ../${target_name}_masked.genome.fasta -bed - -fo ${target_name}_rnd1.all.maker.transcripts.fasta

# remove the duplicated sequence
input_file="BRA1306_rnd1.all.maker.transcripts.fasta"
output_file="BRA1306_rnd1.unique.maker.transcripts.fasta"

# Remove duplicate sequences using awk
awk '/^>/ { if (seq) { print seq; seq=""; } print; next } { seq = seq $0 } END { print seq }' "$input_file" | \
awk '!seen[$0]++' > "$output_file"


# use busco to do evaluation
conda activate busco


PPN=64
target_name=BRA1306
busco -l arthropoda -m trans -i ${target_name}_rnd1.unique.maker.transcripts.fasta --out ${target_name}_maker1_trans -c ${PPN} -f -l /home/lynn/anaconda3/envs/busco/busco_downloads/lineages/arthropoda_odb10 --offline
