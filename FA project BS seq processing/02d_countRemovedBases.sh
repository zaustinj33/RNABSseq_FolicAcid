#!/bin/bash

pushd $2/Code/$1/TecPap

RAW1=$(find /projects/epigenomicslab/RNA_project/ -name $1_1.fq.gz)
RAW2=$(find /projects/epigenomicslab/RNA_project/ -name $1_2.fq.gz)

echo -e "#!/bin/bash
#SBATCH -p normal_q
#SBATCH -A davidxie_lab
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-user zaustinj@vt.edu
#SBATCH --mail-type=END

module reset
module load seqtk

cd $2/raw_data/$1

seqtk fqchk $RAW1 | grep "ALL" > $1_base_counts.txt
seqtk fqchk $RAW2 | grep "ALL" >> $1_base_counts.txt

seqtk fqchk $1_1_trim.fq.gz | grep "ALL" >> $1_base_counts.txt
seqtk fqchk $1_2_trim.fq.gz | grep "ALL" >> $1_base_counts.txt

seqtk fqchk $1_1_qual.fq.gz | grep "ALL" >> $1_base_counts.txt
seqtk fqchk $1_2_qual.fq.gz | grep "ALL" >> $1_base_counts.txt

printf '%s\n' $1 $1 $1 $1 $1 $1 > $1_name_col.txt

paste $1_name_col.txt ../../read_cols.txt > $1_label_table.txt
paste $1_label_table.txt $1_base_counts.txt > $1_base_counts_table.txt

rm $1_base_counts.txt
rm $1_name_col.txt

" >> $1_nucleotideCount.sbatch

sbatch $1_nucleotideCount.sbatch
