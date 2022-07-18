#!/bin/bash

mkdir -p $2/Code/$1/Ch2_code
mkdir -p $2/result/$1
mkdir -p $2/result/$1/processed
mkdir -p $2/raw_data/$1/meRanGhUnaligned
mkdir -p $2/raw_data/$1/meRanTUnaligned

pushd $2/Code/$1/
echo -e "#!/bin/bash
#SBATCH -p normal_q
#SBATCH -A davidxie_lab
#SBATCH --nodes=1 --cpus-per-task=40
#SBATCH --exclusive
#SBATCH --mem=100GB
#SBATCH --time=36:00:00
#SBATCH --mail-user zaustinj@vt.edu
#SBATCH --mail-type=END

# Adapter trim #

module reset
module load fastp

cd $2/working_data/$1

module reset
module load HISAT2
gunzip $1_1_qual.fq.gz
gunzip $1_2_qual.fq.gz

echo 'Start meRanGh'
meRanGh align -o $2/result/$1/ -un -ud $2/raw_data/$1/meRanGhUnaligned -f $1_2_qual.fq -r $1_1_qual.fq -t 40 -fmo -ds -S $2/result/$1/$1_meRanGh_genomeMap.sam -ds -MM -id /projects/epigenomicslab/Annotation/mm10_meRanGh/BSgenomeIDX -GTF /projects/epigenomicslab/Annotation/mm10.ensGene.for.RNABS.gtf
echo 'Finished meRanGh'

gzip $1_1_qual.fq.gz
gzip $1_2_qual.fq.gz

" > $1_R2pos_GNMap.sbatch

sbatch $1_R2pos_GNMap.sbatch
