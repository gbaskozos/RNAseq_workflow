#!/bin/bash

mkdir /data/ndcn-rnaseq/ndcn0569/RNAseq/star_out

LANE1="/data/ndcn-rnaseq/ndcn0569/RNAseq/bsg-ftp.well.ox.ac.uk/180726_K00150_0353_BHWLCMBBXX"

#LANE1
rm -rf $LANE1/files.txt

POS=$((${#LANE1} + 1))
LENGTH=19

for FILE in $LANE1/*.fastq.gz
do
echo ${FILE:POS:LENGTH} >> $LANE1/files.txt
done

sort $LANE1/files.txt | uniq

for FILE in $(sort $LANE1/files.txt | uniq)
do
echo "Mapping reads for (ENCODE parameters)" $FILE
rm -fR /data/ndcn-rnaseq/ndcn0569/RNAseq/STAR_TEMP
$HOME/bin/STAR --runThreadN 15 --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --genomeDir /data/ndcn-rnaseq/ndcn0569/hg_genome_index_75 --readFilesIn "${LANE1}/${FILE}_1.fastq.gz" "${LANE1}/${FILE}_2.fastq.gz"  --readFilesCommand gunzip -c --outFileNamePrefix "$DATA/RNAseq/star_out/$FILE" --outTmpDir $DATA/RNAseq/STAR_TEMP --outSAMtype BAM SortedByCoordinate
echo "Waiting BAMsort, deleting temp files"
rm -fR /data/ndcn-rnaseq/ndcn0569/RNAseq/STAR_TEMP

done

wait

echo "finished mapping"

for file in $STAR_OUT/*Log.final.out
do

awk -F "\t" 'NR == 10 {print FILENAME "\t" $2}' $file | sed 's/%//' >> $DATA/star_out/uniquelly_mapped.txt


done
