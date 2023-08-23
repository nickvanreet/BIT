#!/bin/bash
#SBATCH --job-name=TBRunique            # Job name
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

python3 tsne_loop_CN2.py