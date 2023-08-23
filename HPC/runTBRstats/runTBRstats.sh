#!/bin/bash
#SBATCH --job-name=TBRreco            # Job name
#SBATCH --output=out_%j.log           # Output log file with Job ID
#SBATCH --error=error_%j.log          # Error log file with Job ID
#SBATCH --ntasks=1                    # Run on a single node
#SBATCH --cpus-per-task=28            # Number of CPU cores per task
#SBATCH --mem=200G                    # Total memory limit
#SBATCH --time=24:00:00                # Time limit hrs:min:sec
#SBATCH --mail-type=FAIL              # Type of email notification- FAIL
#SBATCH --mail-user=nvanreet@itg.be   # Your email

# Load Python Module
#module load BioTools-Python/2020a.00-intel-2020a-IntelPython3
module load BioTools-Python/2020a.00-intel-2020a-Python-3.8.3
#module load BioTools

# Python Script Location
PYTHON_SCRIPT="/scratch/antwerpen/206/vsc20643/scripts/HPC/runRECO7/TBRstats.py"

echo $PYTHON_SCRIPT
echo "Found script"

# Define Root Fasta Directory
ROOT_FASTA_DIR="/scratch/antwerpen/206/vsc20643/TBRhunt/FASTA/"

echo $ROOT_FASTA_DIR
echo "Found input"

# Define Root RECO Directory
ROOT_RECO_DIR="/scratch/antwerpen/206/vsc20643/TBRhunt/TBRstats/"

echo $ROOT_RECO_DIR
echo "Found output"

# Execute Python Script
python3 $PYTHON_SCRIPT -d $ROOT_FASTA_DIR -o $ROOT_RECO_DIR

echo "Lets kick some shit"