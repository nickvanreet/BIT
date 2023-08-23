#!/bin/bash
#SBATCH --job-name=TBRloop
#SBATCH --output=out_%j.log
#SBATCH --error=error_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=200G
#SBATCH --time=72:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nvanreet@itg.be

reference_genome=/scratch/antwerpen/206/vsc20643/trypanosoma_genomes/data/refgenome/TriTrypDB-63_TbruceiTREU927_Genome_withTBR3.fasta
threads=28

module load BioTools
module load BWA

root_bwa_dir="/scratch/antwerpen/206/vsc20643/TBRhunt/BWA/"
root_fastq_dir="/data/antwerpen/grp/aitg/pbuscher/Trypanosoma/rawreads/BATCH3"

# Loop through subdirectories (strains)
for strain_dir in $root_fastq_dir/*/ ; do
  strain=$(basename $strain_dir)
  echo "Processing strain: $strain"
  
  # Loop through fastq pairs within each strain
  for fastq_file_1 in $strain_dir/*_1.fq.gz ; do
    fastq_file_2=${fastq_file_1/_1.fq.gz/_2.fq.gz}
    
    bwa_dir=${root_bwa_dir}/${strain}
    mkdir -p ${bwa_dir}

    bam_file_prefix=${bwa_dir}/${strain}

    bwa mem \
    -R "@RG\tID:${strain}\tSM:${strain}\tPL:ILLUMINA" \
    -t $threads \
    $reference_genome \
    $fastq_file_1 \
    $fastq_file_2 \
    | samtools sort -@ $threads -o ${bam_file_prefix}.bam

    samtools index ${bam_file_prefix}.bam
    samtools flagstat -@ ${threads} ${bam_file_prefix}.bam > ${bam_file_prefix}.flagstat

    samtools view -b ${bam_file_prefix}.bam TBR3 > ${bam_file_prefix}.TBR.bam
    samtools index ${bam_file_prefix}.TBR.bam
    samtools flagstat -@ ${threads} ${bam_file_prefix}.TBR.bam > ${bam_file_prefix}.TBR.flagstat

    rm ${bam_file_prefix}.bam
  done
done

echo "Finished BWA"
