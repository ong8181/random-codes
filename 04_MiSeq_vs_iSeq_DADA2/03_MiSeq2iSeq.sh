####
#### Check how iSeq Q-score influences DADA2 results
#### 2020.7.22 Ushio
####

FASTQ_FOLDER1="02_DemultiplexedOut_MiSeqStyle"
FASTQ_FOLDER2="02_DemultiplexedOut_iSeqStyle"

# 1. Should prepare MiSeq-style fastq files
mkdir $FASTQ_FOLDER2
cp $FASTQ_FOLDER1/*.gz $FASTQ_FOLDER2

# 2. Load function to convert MiSeq-style fastq to iSeq-style
source 00_Helpers/MiSeq2iSeq.sh

# 3. Convert MiSeq-Style fastq files to iSeq-Style
cd $FASTQ_FOLDER2
pigz -d *.gz
for file in *; do miseq2iseq $file; done
pigz *.fastq
cd ..



