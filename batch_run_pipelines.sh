NCORE=24
MEM=100

samples=`ls $DATA_PATH | grep -oP '.+(?=_Transcriptome)'`
ref=`readlink -f refdata-gex-GRCh38-2020-A`
for sample in $samples; do
    bsub -q verylong -n $NCORE -R "rusage[mem=${MEM}G]" -J cellranger_$sample scripts/run_cellranger.sh outputs $sample sample_id_barcode.txt $NCORE $MEM $ref
done