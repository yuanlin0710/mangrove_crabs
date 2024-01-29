fasta_path=/home/lynn/data/fyp/ABEL_CRAB_DATA
target_name=BRA1313
species_name=Episesarma_versicolor
PPN=64
/home/lynn/program/PASApipeline-master/Launch_PASA_pipeline.pl -c pasa.alignAssembly.BRA1313_Mar15update.txt -C -R -g ${target_name}_masked.merge.fasta --TDN ${target_name}-all_tissues.Trinity.uniq.raname.txt --ALIGNERS minimap2 -t ${target_name}_genomic_masked_mapped+trinity_transcripts.fa --CPU $PPN
