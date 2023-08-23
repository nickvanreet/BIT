#!/bin/bash
#SBATCH --job-name=runFLAGSTAT
#SBATCH --output=out_%j.log
#SBATCH --error=error_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=200G
#SBATCH --time=1:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nvanreet@itg.be

# Output file
out_file="/scratch/antwerpen/206/vsc20643/TBRhunt/FLAGSTAT/combined_flagstat_report.csv"
FLAGSTAT_dir="/scratch/antwerpen/206/vsc20643/TBRhunt/BWA/"

# Header for the csv file
echo "Strain,Type,Total,QC_passed,QC_failed,Secondary,Supplementary,Duplicates,Mapped,Mapped_percent,Paired,Read1,Read2,Properly_paired,Properly_paired_percent,Self_and_mate_mapped,Singletons,Singletons_percent,Mate_mapped_diff_chr,Mate_mapped_diff_chr_mapQ_ge_5" > $out_file

# Loop over subdirectories
for dir in ${FLAGSTAT_dir}*/; do
    # Get strain name from directory name
    strain=$(basename $dir)

    # Parse flagstat file for this strain
    for flagstat_file in ${dir}*.flagstat; do
        # If the .flagstat file doesn't exist, skip this iteration
        if [ ! -f "$flagstat_file" ]; then
            continue
        fi

        # Determine type based on file name
        if [[ $flagstat_file == *".TBR.flagstat" ]]; then
            type="TBR"
        else
            type="Normal"
        fi

        # Use awk to parse flagstat file
        parsed_data=$(awk '
            /in total/ {total=$1; qc_passed=$1; qc_failed=$3}
            /secondary/ {secondary=$1}
            /supplementary/ {supplementary=$1}
            /duplicates/ {duplicates=$1}
            /mapped/ {mapped=$1; mapped_percent=$3}
            /paired in sequencing/ {paired=$1}
            /read1/ {read1=$1}
            /read2/ {read2=$1}
            /properly paired/ {properly_paired=$1; properly_paired_percent=$3}
            /with itself and mate mapped/ {self_and_mate_mapped=$1}
            /singletons/ {singletons=$1; singletons_percent=$3}
            /with mate mapped to a different chr/ {mate_mapped_diff_chr=$1}
            /with mate mapped to a different chr \(mapQ>=5\)/ {mate_mapped_diff_chr_mapQ_ge_5=$1}
            END {print total, qc_passed, qc_failed, secondary, supplementary, duplicates, mapped, mapped_percent, paired, read1, read2, properly_paired, properly_paired_percent, self_and_mate_mapped, singletons, singletons_percent, mate_mapped_diff_chr, mate_mapped_diff_chr_mapQ_ge_5}
        ' $flagstat_file)

        # Add strain name and type to the start of the parsed data
        csv_row="$strain,$type,$parsed_data"

        # Append row to csv file
        echo $csv_row >> $out_file
    done
done
