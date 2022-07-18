#!/bin/bash

#SBATCH -p normal_q
#SBATCH -A davidxie_lab
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --time=1:00:00
# SBATCH --mail-user zaustinj@vt.edu
# SBATCH --mail-type=END

# How to Submit: sbatch [this file name] $PWD $file1 $file2

#OGcode: 
#perl /groups/ECBL/project_data/VTCRI_Methylome/DMA/DMR_scripts/Multi_Thread_DMS_Processor.pl -r /apps/easybuild/software/pegasus-sandy_bridge/R/3.5.0-foss-2018a-X11-20180131/bin/ /home/bsharmi6/RNA_methylation/G3lowFARNABS_bsharmi.txt /home/bsharmi6/RNA_methylation/G1controlRNABS_bsharmi.txt  -o /home/bsharmi6/RNA_methylation/DMS_output/ -t /home/bsharmi6/RNA_methylation/DMS_temp/

# Add any modules you might require:
module load R/4.0.2-foss-2020a
mkdir -p $1"/DMS/DMS_output"

perl /projects/epigenomicslab/zaustinj/RNA_BSnew/miniScripts/Multi_Thread_DMS_Processor.pl \
-r /apps/easybuild/software/tinkercliffs-rome/R/4.0.2-foss-2020a/bin "/beegfs/projects/epigenomicslab/RNA_project/mouse_NSC_FA_project/FA_BS/DMS/"$2 "/beegfs/projects/epigenomicslab/RNA_project/mouse_NSC_FA_project/FA_BS/DMS/"$3 \
-o $1"/DMS/DMS_output/tLF_x_tMF" \
-t $1"/DMS/DMS_temp/"

directory=$1"/DMS/DMS_output/tLF_x_tMF"
outfile="$(ls -1 $directory)"


Rscript /projects/epigenomicslab/project_data/VTCRI_Methylome/DMA/DMR_scripts/DMS.p.adj.bsharmi.r $directory/$outfile RNA_meth 10

