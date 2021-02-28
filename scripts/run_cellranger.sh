#!/bin/bash

ANA_PATH=$1
SAMPLE=$2
SAMPLE_ID_BARCODE_FILE=$3
NCORE=$4
MEM=$5
REF=$6

mkdir -p $ANA_PATH/$SAMPLE
echo "id,name,read,pattern,sequence,feature_type" > $ANA_PATH/$SAMPLE/feature_ref.csv
cat $SAMPLE_ID_BARCODE_FILE | grep $SAMPLE | cut -f 2,3 | awk '{split($0,a,"\t"); print a[1]","a[1]",R2,^(BC),"a[2]",Antibody Capture"}' >> $ANA_PATH/$SAMPLE/feature_ref.csv
echo "fastqs,sample,library_type" > $ANA_PATH/$SAMPLE/lib.csv
fastq_path=$DATA_PATH/${SAMPLE}_Antibody_Multiplex_CITE
fastq_head=`ls $fastq_path | grep -oP '.+(?=_S[0-9]+_L[0-9]+_R1_001.fastq.gz)'`
echo "$fastq_path,$fastq_head,Antibody Capture" >> $ANA_PATH/$SAMPLE/lib.csv
fastq_path=$DATA_PATH/${SAMPLE}_Transcriptome
fastq_head=`ls $fastq_path | grep -oP '.+(?=_S[0-9]+_L[0-9]+_R1_001.fastq.gz)'`
echo "$fastq_path,$fastq_head,Gene Expression" >> $ANA_PATH/$SAMPLE/lib.csv
cd $ANA_PATH/$SAMPLE && cellranger count --id=$SAMPLE --transcriptome=$REF --feature-ref=feature_ref.csv --libraries=lib.csv --localcores=$NCORE --localmem=$MEM --expect-cells=10000 --chemistry=SC3Pv3 > $ANA_PATH/$SAMPLE/log.txt 2>&1