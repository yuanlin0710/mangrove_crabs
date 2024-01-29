target_name=bra1313
TransDecoder.LongOrfs -t ${target_name}_sqlite_output.assemblies.fasta
TransDecoder.Predict -t ${target_name}_sqlite_output.assemblies.fasta
cdna_alignment_orf_to_genome_orf.pl ${target_name}_sqlite_output.assemblies.fasta.transdecoder.gff3 ${target_name}_sqlite_output.pasa_assemblies.gff3 ${target_name}_sqlite_output.assemblies.fasta > ${target_name}_sqlite_output.assemblies.fasta.transdecoder.genome.gff3
