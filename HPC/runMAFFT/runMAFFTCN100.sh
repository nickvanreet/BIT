#!/bin/bash
#SBATCH --job-name=runMAFFT            # Job name
#SBATCH --output=out_%j.log           # Output log file with Job ID
#SBATCH --error=error_%j.log          # Error log file with Job ID
#SBATCH --ntasks=1                    # Run on a single node
#SBATCH --cpus-per-task=28            # Number of CPU cores per task
#SBATCH --mem=200G                    # Total memory limit
#SBATCH --time=24:00:00                # Time limit hrs:min:sec
#SBATCH --mail-type=FAIL              # Type of email notification- FAIL
#SBATCH --mail-user=nvanreet@itg.be   # My email

threads=28

# Load necessary module for MAFFT
module load BioTools-Python/2020a.00-intel-2020a-IntelPython3
module load MAFFT/7.471-intel-2020a-with-extensions

# Directory for MAFFT input and output
MAFFT_input="//scratch/antwerpen/206/vsc20643/TBRhunt/UNIQUE/CN_100/unique_sequences_CN100.fasta"
MAFFT_output="/scratch/antwerpen/206/vsc20643/TBRhunt/MAFFT/"
MAFFT_output_aln="${MAFFT_output}/aligned_unique_sequences_CN100.fasta"

# Run MAFFT on the specified FASTA files
./simpleMAFFT_CN.py -i "$MAFFT_input" -o "$MAFFT_output_aln" -t "$threads"
echo "Finished MAFFT alignment"