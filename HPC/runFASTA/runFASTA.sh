#!/bin/bash
#SBATCH --job-name=runFASTA
#SBATCH --output=out_%j.log
#SBATCH --error=error_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nvanreet@itg.be

# Extract the FASTA files from TBR.bam
module load BioTools-Python/2020a.00-intel-2020a-IntelPython3

# Directories for process_bam.py
root_bwa_dir="/scratch/antwerpen/206/vsc20643/TBRhunt/BWA/"
root_fasta_dir="/scratch/antwerpen/206/vsc20643/TBRhunt/FASTA/"

# Loop through subdirectories (strains)
for strain_dir in $root_bwa_dir/*/ ; do
  strain=$(basename $strain_dir)
  echo "Processing strain: $strain"

  # Loop through BWA results within each strain
  fasta_dir=${root_fasta_dir}/${strain}
  # Create directory if it doesn't exist
  mkdir -p ${fasta_dir}

  # Python script to process BAM files
  ./process_bam_HPC.py $strain_dir $fasta_dir

  echo "Finished processing BAM files to FASTA"
done
