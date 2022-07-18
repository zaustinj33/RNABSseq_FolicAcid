#!/bin/bash

mkdir -p $2/Code/$1/
mkdir -p $2/working_data/$1
mkdir -p $2/result/$1

pushd $2/Code/$1/
RAW1=$(find /projects/epigenomicslab/RNA_project/ -name $1_1.fq.gz)
RAW2=$(find /projects/epigenomicslab/RNA_project/ -name $1_2.fq.gz)

echo -e "#!/bin/bash
#SBATCH -p normal_q
#SBATCH -A davidxie_lab
#SBATCH --nodes=1 --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --mail-user zaustinj@vt.edu
#SBATCH --mail-type=END

cd $2/raw_data/$1

module reset
module load fastp

cd $2/working_data/$1

# Too short or no mate, remove adapters, remove polyx tails
fastp -w 16 -Q -l 50 --trim_poly_x --poly_x_min_len 10 -i $2/raw_data/$1/$1_1.fq.gz -I $2/raw_data/$1/$1_2.fq.gz --out1 $1_1_trim.fq.gz --out2 $1_2_trim.fq.gz --failed_out $1_failed_length.fq.gz -j $1_failed_length.json -h $1_failed_length.html --detect_adapter_for_pe --overlap_diff_percent_limit 25

# Remove low quality bases; window trim
fastp -w 16 -q 25 -5 -3 -M 25 -i $1_1_trim.fq.gz -I $1_2_trim.fq.gz --out1 $1_1_qual.fq.gz --out2 $1_2_qual.fq.gz --failed_out $1_failed_qual.fq.gz -j $1_failed_qual.json -h $1_failed_qual.html --detect_adapter_for_pe


" > $1_sequenceTrim.sbatch

sbatch $1_sequenceTrim.sbatch
