#!/bin/bash

mkdir -p $2/Code/$1/Ch4_code
mkdir -p $2/CallResult/$1/
pushd $2/Code/$1/Ch4_code

GNconvrate=$(cat $2/Code/$1/Ch3_code/$1_GN_conv_rate.txt)

echo -e "#!/bin/bash
#SBATCH -p normal_q
#SBATCH -A davidxie_lab
#SBATCH --nodes=1 --cpus-per-task=40
#SBATCH --exclusive
#SBATCH --time=100:00:00
#SBATCH --mail-user zaustinj@vt.edu
#SBATCH --mail-type=END

mkdir -p \$TMPDIR/ZJ_tmp


## Genome
# Call sites
#meRanCall -p 40 -fs5 6 -fs3 6 -rs5 6 -rs3 6 -o \$TMPDIR/ZJ_tmp/$1_Genome10xCall.txt -bam $2/result/$1/$1_meRanGh_genomeMapMerged_sorted.bam -f /projects/epigenomicslab/Annotation/mm10.for.RNABS.fa -mBQ 30 -gref -rl 150 -sc 10 -cr 1 -mr 0.00001 -mcov 10
#meRanCall -p 40 -fs5 6 -fs3 6 -rs5 6 -rs3 6 -o \$TMPDIR/ZJ_tmp/$1_Genome10xCall_3_Cutoff.txt -bam $2/result/$1/$1_3_Ccutoff_PE.bam -f /projects/epigenomicslab/Annotation/mm10.for.RNABS.fa -mBQ 30 -gref -rl 150 -sc 10 -cr 1 -fdr 0.01 -mr 0.00001 -mcov 10
#meRanCall -p 40 -fs5 6 -fs3 6 -rs5 6 -rs3 6 -o \$TMPDIR/ZJ_tmp/$1_Genome10xCall_3_Cutoff_FDR.txt -bam $2/result/$1/$1_3_Ccutoff_PE.bam -f /projects/epigenomicslab/Annotation/mm10.for.RNABS.fa -mBQ 30 -gref -rl 150 -sc 10 -cr $GNconvrate -fdr 0.05 -mr 0.00001 -mcov 10

#mv \$TMPDIR/ZJ_tmp/$1_Genome10xCall.txt $2/CallResult/$1/$1_Genome10xCall.txt
#mv \$TMPDIR/ZJ_tmp/$1_Genome10xCall_3_Cutoff.txt $2/CallResult/$1/$1_Genome10xCall_3_Cutoff.txt
#mv \$TMPDIR/ZJ_tmp/$1_Genome10xCall_3_Cutoff_FDR_FDR_0.05.txt $2/CallResult/$1/$1_Genome10xCall_3_Cutoff_FDR_FDR_0.05.txt

# Annotate sites
#sed -i 's/chrM/chrMT/g' $2/CallResult/$1/$1_Genome10xCall.txt
#sed -i 's/chr//' $2/CallResult/$1/$1_Genome10xCall.txt

#sed -i 's/chrM/chrMT/g' $2/CallResult/$1/$1_Genome10xCall_3_Cutoff.txt
#sed -i 's/chr//' $2/CallResult/$1/$1_Genome10xCall_3_Cutoff.txt 

#meRanAnnotate -t $2/CallResult/$1/$1_Genome10xCall.txt -ensGTF -g /projects/epigenomicslab/Annotation/Mus_musculus.GRCm38.96.gtf -o $2/CallResult/$1/$1_Genome10xCall_annotate.txt -f 'gene'
#meRanAnnotate -t $2/CallResult/$1/$1_Genome10xCall_3_Cutoff.txt -ensGTF -g /projects/epigenomicslab/Annotation/Mus_musculus.GRCm38.96.gtf -o $2/CallResult/$1/$1_Genome10xCall_3_Cutoff_annotate.txt -f 'gene'

# Only depth sites
meRanCall -p 40 -o \$TMPDIR/ZJ_tmp/$1_Genome10xCall_3_Cutoff_0ML.txt -bam $2/result/$1/$1_3_Ccutoff_PE.bam -f /projects/epigenomicslab/Annotation/mm10.for.RNABS.fa -mBQ 30 -gref -rl 150 -sc 10 -cr 1 -mr 0 -mcov 10 -regions $2/All_m5C.bed
mv \$TMPDIR/ZJ_tmp/$1_Genome10xCall_3_Cutoff_0ML.txt $2/CallResult/$1/$1_Genome10xCall_3_Cutoff_0ML.txt



" > $1_meRanCallCutoff_0ML.sbatch

sbatch $1_meRanCallCutoff_0ML.sbatch
