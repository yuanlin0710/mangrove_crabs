MODE=coassembly
P_NAME=EV
S_FILE=/home/lynn/data/fyp_gut_metagenome/squeezemeta/sample_EV
FQ_DIR=/home/lynn/data/fyp_gut_metagenome/trimmed_read
EXTDB_DIR=/home/lynn/data/fyp_gut_metagenome/squeezemeta/mydb.list
THREAD_NUM=64
EXT_ASSEMBLY=/home/lynn/data/fyp_gut_metagenome/assembly_external_Tom/01.EV1_3.fasta

SqueezeMeta.pl -m $MODE -p $P_NAME -s $S_FILE -f $FQ_DIR --extdb $EXTDB_DIR -t $THREAD_NUM -extassembly $EXT_ASSEMBLY





(SqueezeMeta) lynn@poweredge:~/data/fyp_gut_metagenome/squeezemeta$ more sample_EV
EV1     EV1.trim.R1.fq.gz       pair1
EV1     EV1.trim.R2.fq.gz       pair2
EV2     EV2.trim.R1.fq.gz       pair1
EV2     EV2.trim.R2.fq.gz       pair2
EV3     EV3.trim.R1.fq.gz       pair1
EV3     EV3.trim.R2.fq.gz       pair2
