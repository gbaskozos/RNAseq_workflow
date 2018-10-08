#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=8:00:00
#SBATCH --job-name=count_CTS
#SBATCH --mail-type=ALL
#SBATCH --mail-user=georgios.baskozos@ndcn.ox.ac.uk

python="/home/ndcn-rnaseq/ndcn0569/python/bin/python"

BAM_FILES="/data/ndcn-rnaseq/ndcn0569/RNAseq/star_out/merged_bams"


for file in $BAM_FILES/*1.sorted.bam
do
echo $file

srun $python -m HTSeq.scripts.count -f bam -r name -s reverse -m intersection-nonempty $file /data/ndcn-rnaseq/ndcn0569/hg_genome/Homo_sapiens.GRCh38.88.chr.gtf > $file.counts.tab           
        
done &

for file in $BAM_FILES/*2.sorted.bam
do
echo $file

srun $python -m HTSeq.scripts.count -f bam -r name -s reverse -m intersection-nonempty $file /data/ndcn-rnaseq/ndcn0569/hg_genome/Homo_sapiens.GRCh38.88.chr.gtf > $file.counts.tab           
        
done &

for file in $BAM_FILES/*3.sorted.bam
do
echo $file

srun $python -m HTSeq.scripts.count -f bam -r name -s reverse -m intersection-nonempty $file /data/ndcn-rnaseq/ndcn0569/hg_genome/Homo_sapiens.GRCh38.88.chr.gtf > $file.counts.tab           
        
done &

for file in $BAM_FILES/*4.sorted.bam
do
echo $file

srun $python -m HTSeq.scripts.count -f bam -r name -s reverse -m intersection-nonempty $file /data/ndcn-rnaseq/ndcn0569/hg_genome/Homo_sapiens.GRCh38.88.chr.gtf > $file.counts.tab           
        
done &
   
       
for file in $BAM_FILES/*5.sorted.bam
do
echo $file

srun $python -m HTSeq.scripts.count -f bam -r name -s reverse -m intersection-nonempty $file /data/ndcn-rnaseq/ndcn0569/hg_genome/Homo_sapiens.GRCh38.88.chr.gtf > $file.counts.tab           
        
done &


for file in $BAM_FILES/*6.sorted.bam
do
echo $file

srun $python -m HTSeq.scripts.count -f bam -r name -s reverse -m intersection-nonempty $file /data/ndcn-rnaseq/ndcn0569/hg_genome/Homo_sapiens.GRCh38.88.chr.gtf > $file.counts.tab           
        
done &


for file in $BAM_FILES/*7.sorted.bam
do
echo $file

srun $python -m HTSeq.scripts.count -f bam -r name -s reverse -m intersection-nonempty $file /data/ndcn-rnaseq/ndcn0569/hg_genome/Homo_sapiens.GRCh38.88.chr.gtf > $file.counts.tab           
        
done &


for file in $BAM_FILES/*8.sorted.bam
do
echo $file

srun $python -m HTSeq.scripts.count -f bam -r name -s reverse -m intersection-nonempty $file /data/ndcn-rnaseq/ndcn0569/hg_genome/Homo_sapiens.GRCh38.88.chr.gtf > $file.counts.tab           
        
done &


for file in $BAM_FILES/*9.sorted.bam
do
echo $file

srun $python -m HTSeq.scripts.count -f bam -r name -s reverse -m intersection-nonempty $file /data/ndcn-rnaseq/ndcn0569/hg_genome/Homo_sapiens.GRCh38.88.chr.gtf > $file.counts.tab           
        
done &


for file in $BAM_FILES/*0.sorted.bam
do
echo $file

srun $python -m HTSeq.scripts.count -f bam -r name -s reverse -m intersection-nonempty $file /data/ndcn-rnaseq/ndcn0569/hg_genome/Homo_sapiens.GRCh38.88.chr.gtf > $file.counts.tab           
        
done &

wait


