#!/bin/bash

#SBATCH --ntasks=1 --cpus-per-task=28
#SBATCH --time=20:00:00
# #SBATCH --job-name=bwa_run

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nvanreet@itg.be

#SBATCH -o stdout.%j
#SBATCH -e stderr.%j


## load necessary module ##
module load BioTools
module load BWA

#### bwa index /user/antwerpen/206/vsc20643/scratch/trypanosoma_genomes/data/refgenome/TriTrypDB-54_TbruceiTREU927_Genome_withTBR.fasta
#### zeker niet in script laten , eeerst in terminal intikken zonder SBATCH

## set parameters
# reference_genome=/user/antwerpen/205/vsc20587/scratch/trypanosoma_genomes/data/refgenome/TriTrypDB-63_TbruceiTREU927_Genome.fasta
reference_genome=/user/antwerpen/206/vsc20643/scratch/trypanosoma_genomes/data/refgenome/TriTrypDB-54_TbruceiTREU927_Genome_withTBR.fasta
bwa_dir=/user/antwerpen/206/vsc20643/scratch/trypanosoma_genomes/results/bwa/
fastq_dir=/user/antwerpen/206/vsc20643/scratch/trypanosoma_genomes/results/fastq/
threads=28


## fastq_file_R1: should be read in via the parameter setting but below
## an example fastq_file_R1 file
fastq_file_1=/data/antwerpen/grp/aitg/pbuscher/Trypanosoma/rawreads/BATCH1/348BT/DP8400009737BL_L01_SP2004030342_1.fq.gz
fastq_file_2=/data/antwerpen/grp/aitg/pbuscher/Trypanosoma/rawreads/BATCH1/348BT/DP8400009737BL_L01_SP2004030342_2.fq.gz
# fastq_file_1=/user/antwerpen/205/vsc20587/aitg_data/pbuscher/Trypanosoma/rawreads/BATCH4/BIP42/BIP42_1.fq.gz
# fastq_file_2=/user/antwerpen/205/vsc20587/aitg_data/pbuscher/Trypanosoma/rawreads/BATCH4/BIP42/BIP42_2.fq.gz
fastq_file_1=${fastq_file_1}
fastq_file_2=${fastq_file_2}

echo "fastq_file_1 & 2..."
echo $fastq_file_1
echo $fastq_file_2

## extract the prefix file and create name 
## for the concatenated fastq file. Extract the 
## sample name
file_prefix_full=${fastq_file_1%_1.fq.gz}
file_prefix_full=${file_prefix_full%_1.fastq.gz}
sample=${file_prefix_full##*/}
fastq_file_concat=${fastq_dir}/${sample}.fastq.gz

## do concatenation
cat ${fastq_file_1} ${fastq_file_2} > ${fastq_file_concat}

## create the bam_file name
bam_file_prefix=${bwa_dir}/${sample}

echo $fastq_file_1
echo $fastq_file_2
echo $sample
echo $bam_file_prefix


## align the reads first to the human genome to filter out the
## native human reads
bwa mem \
	-R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" \
	-t $threads \
	$reference_genome \
	$fastq_file_concat \
	| samtools sort -@ $threads -o ${bam_file_prefix}.bam

## run flagstat for statistics
samtools index ${bam_file_prefix}.bam
samtools flagstat -@ ${threads} ${bam_file_prefix}.bam > ${bam_file_prefix}.flagstat

## extract only the TBR region
samtools view -b ${bam_file_prefix}.bam TBR > ${bam_file_prefix}.TBR.bam
samtools index ${bam_file_prefix}.TBR.bam
samtools flagstat -@ ${threads} ${bam_file_prefix}.TBR.bam > ${bam_file_prefix}.TBR.flagstat

## clean up 
rm ${bam_file_prefix}.bam
rm ${fastq_file_concat}

## run script per fastq file\
# cd /scratch/antwerpen/grp/aitg/arosanas/pvivax_genomes/data/fastq_new
# sbatch --export=fastq_file_1=/scratch/antwerpen/grp/aitg/arosanas/pvivax_genomes/data/fastq_new/ERR5740834_1.fastq.gz /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/BWA_run.sh

# while IFS=, read -r fastq_file_1 fastq_file_2; do 
# 	echo "FastQ_file1: $fastq_file_1 :: FastQ_file2: $fastq_file_2"; 
# 	sbatch --export=fastq_file_1=${fastq_file_1},fastq_file_2=${fastq_file_2} /user/antwerpen/205/vsc20587/scratch/trypanosoma_genomes/bin/BWA_run_TBR.sh
# done < /user/antwerpen/205/vsc20587/scratch/trypanosoma_genomes/results/fastq/fastq_files.csv

## run again for those samples that crashed due to 
## write protected folders in the ait_data of Manon (BATCH 3 & 4)
# while IFS=, read -r fastq_file_1 fastq_file_2; do 
# 	echo "FastQ_file1: $fastq_file_1 :: FastQ_file2: $fastq_file_2"; 
# 	sbatch --export=fastq_file_1=${fastq_file_1},fastq_file_2=${fastq_file_2} /user/antwerpen/205/vsc20587/scratch/trypanosoma_genomes/bin/BWA_run_TBR.sh
# done < /user/antwerpen/205/vsc20587/scratch/trypanosoma_genomes/results/fastq/fastq_files.run2.csv