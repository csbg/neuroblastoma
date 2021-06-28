#!/bin/bash
echo "Processing sample $1"

OMICSTMP="/netscratch/omicstmp/WES"
OMICSTMP_SIF="/media/omicstmp"
RAWDATADIR="/usr/local/AGFORTELNY/DATA_Neuroblastoma/data_raw/COUNT/$1_transcriptome"


echo "Copying files to $OMICSTMP/input_$1"

mkdir -p $OMICSTMP/input_$1
mkdir -p $OMICSTMP/output_$1

cp $RAWDATADIR/possorted_genome_bam.bam $OMICSTMP/input_$1/
cp $RAWDATADIR/possorted_genome_bam.bam.bai $OMICSTMP/input_$1/
cp $RAWDATADIR/filtered_feature_bc_matrix/barcodes.tsv.gz $OMICSTMP/input_$1/
gunzip $OMICSTMP/input_$1/barcodes.tsv.gz


echo "Files copied, starting velocyto"

singularity run \
  --bind $OMICSTMP:$OMICSTMP_SIF \
  scvelo.sif \
  velocyto run \
    -b $OMICSTMP_SIF/input_$1/barcodes.tsv \
    -o $OMICSTMP_SIF/output_$1 \
    -e possorted_genome_bam \
    -m $OMICSTMP_SIF/metadata/hg38_rmsk.gtf \
    -@ 1 \
    -vv \
    $OMICSTMP_SIF/input_$1/possorted_genome_bam.bam \
    $OMICSTMP_SIF/metadata/genes.gtf


echo "Removing input files"

rm -rf $OMICSTMP/input_$1
