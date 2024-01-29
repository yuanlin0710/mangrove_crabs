target_name=BRA1306
gff3_merge -s -d ${target_name}_train1_master_datastore_index.log > ${target_name}.all.maker.gff
fasta_merge -d ${target_name}_train1_master_datastore_index.log
gff3_merge -n -s -d ${target_name}_train1_master_datastore_index.log > ${target_name}.all.maker.noseq.gff
