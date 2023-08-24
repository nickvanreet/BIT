#!/bin/bash
#SBATCH --job-name=SNPunique          # Job name
#SBATCH --output=out_%j.log           # Output log file with Job ID
#SBATCH --error=error_%j.log          # Error log file with Job ID
#SBATCH --ntasks=1                    # Run on a single node
#SBATCH --cpus-per-task=28            # Number of CPU cores per task
#SBATCH --mem=200G                    # Total memory limit
#SBATCH --time=24:00:00               # Time limit hrs:min:sec
#SBATCH --mail-type=FAIL              # Type of email notification- FAIL
#SBATCH --mail-user=nvanreet@itg.be   # Your email

threads=28

# Load necessary module for MAFFT
module load BioTools-Python/2020a.00-intel-2020a-IntelPython3

# File paths
alignment_file="/scratch/antwerpen/206/vsc20643/TBRhunt/MAFFT/aligned_unique_sequences_CN100.fasta"
lookup_file="/scratch/antwerpen/206/vsc20643/TBRhunt/UNIQUE/CN_100/copy_number_lookup_CN100.csv"
strain_abb_file="/scratch/antwerpen/206/vsc20643/TBRhunt/META/strain_subspecies_ab_meta.csv"
output_file="/scratch/antwerpen/206/vsc20643/TBRhunt/SNP/CN_100/SNP_CN100.csv"
aggregate_file="/scratch/antwerpen/206/vsc20643/TBRhunt/SNP/CN_100/SNP_CN100_agg.csv"
strain_file="/scratch/antwerpen/206/vsc20643/TBRhunt/SNP/CN_100/SNP_CN100_strain.csv"

# Run the Python script
python3 SNPunique5iter.py -i $alignment_file -o $output_file -l $lookup_file -s $strain_abb_file -a $aggregate_file -sp $strain_file
