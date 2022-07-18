#!/bin/bash

pushd $2/Code/$1/TecPap

echo -e "#!/bin/bash
#SBATCH -p normal_q
#SBATCH -A davidxie_lab
#SBATCH --nodes=1 --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --mail-user zaustinj@vt.edu
#SBATCH --mail-type=END

## get mbias data
# map files
cd $2/raw_data/$1

module reset
module load Bismark
#bismark -p 4 --non_directional --genome /projects/epigenomicslab/Annotation/mm10_bismark/ --se $1_1_qualFormatted.fq.gz -o $2/mapped_files/$1
#bismark -p 4 --non_directional --genome /projects/epigenomicslab/Annotation/mm10_bismark/ --se $1_2_qualFormatted.fq.gz -o $2/mapped_files/$1
cd $2/mapped_files/$1

module load SAMtools
samtools sort -@ 8 $1_1_qualFormatted_bismark_bt2.bam > $1_1_SEsort.bam
samtools sort -@ 8 $1_2_qualFormatted_bismark_bt2.bam > $1_2_SEsort.bam
samtools index $1_1_SEsort.bam
samtools index $1_2_SEsort.bam

#module load SAMtools
#samtools sort -n -@ 8 $1_deduplicated.bam > $1_deduplicated_sort.bam
#samtools index $1_deduplicated_sort.bam
#module reset
#module load Bismark
bismark_methylation_extractor -s --mbias_only --parallel 10 --genome_folder /projects/epigenomicslab/Annotation/mm10_bismark/ $1_1_SEsort.bam
bismark_methylation_extractor -s --mbias_only --parallel 10 --genome_folder /projects/epigenomicslab/Annotation/mm10_bismark/ $1_2_SEsort.bam


" >> $1_mBias_SE_call.sbatch

sbatch $1_mBias_SE_call.sbatch
