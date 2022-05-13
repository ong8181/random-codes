####
#### Demultiplexing, Primer trimming and other preprocessing
#### Re-analyze O'Donnell et al. (2016)
#### 2022.05.12
####

#---------------------------------------------------------------------------------------#
# REQIRED: fastp (https://github.com/OpenGene/fastp)
# REQIRED: seqkit (https://bioinf.shenwei.me/seqkit)

# -------------------------------------------- #
# Step 0. Demultiplex and preparations
# -------------------------------------------- #
# Set parameters
WORKING_DIR=/Users/ushio/Desktop/ODonnell_2016/Analysis
TAG_DATA1=sampledata/tag_info1.csv
TAG_DATA2=sampledata/tag_info2.csv
TAG_DATA3=sampledata/tag_info3.csv
TAG_DATA4=sampledata/tag_info4.csv
TAG_DATA5=sampledata/tag_info5.csv
DEMILTIPLEX_DIR=seqdata_demultiplexed

#SEQ_MODE=PE_DualID
SEQ_DIR=seqdata
SEQ_FILE1_R1=library1_S1_L001_R1_001.fastq.gz
SEQ_FILE1_R2=library1_S1_L001_R2_001.fastq.gz
SEQ_FILE2_R1=library2_S2_L001_R1_001.fastq.gz
SEQ_FILE2_R2=library2_S2_L001_R2_001.fastq.gz
SEQ_FILE3_R1=RK16sTS1_S1_L001_R1_001.fastq.gz
SEQ_FILE3_R2=RK16sTS1_S1_L001_R2_001.fastq.gz
SEQ_FILE4_R1=RK16sTS2_S2_L001_R1_001.fastq.gz
SEQ_FILE4_R2=RK16sTS2_S2_L001_R2_001.fastq.gz
SEQ_FILE5_R1=RK16sTS3_S3_L001_R1_001.fastq.gz
SEQ_FILE5_R2=RK16sTS3_S3_L001_R2_001.fastq.gz

# ------------------------------------------------------------------------------------- #
# Demultiplexing using cutadapt
# ------------------------------------------------------------------------------------- #
# Move to the working directory
cd ${WORKING_DIR}


# ------------------------------------ #
# Generate tag fastq files
# ------------------------------------ #
## ----------------------- For sequence data 1
count=0
while read row || [ -n "${row}" ]; do
  if ((count > 0)); then
    sample_name=$(echo ${row} | cut -d , -f 1)
    tag1=$(echo ${row} | cut -d , -f 2)
    tag2=$(echo ${row} | cut -d , -f 3)
    echo -e ">"${sample_name}"\n""^"${tag1} >> sampledata/seq1_R1_index.fasta
    echo -e ">"${sample_name}"\n""^"${tag2} >> sampledata/seq1_R2_index.fasta
  fi
  count=`expr $count + 1`
done < ${TAG_DATA1}
## ----------------------- For sequence data 2
count=0
while read row || [ -n "${row}" ]; do
  if ((count > 0)); then
    sample_name=$(echo ${row} | cut -d , -f 1)
    tag1=$(echo ${row} | cut -d , -f 2)
    tag2=$(echo ${row} | cut -d , -f 3)
    echo -e ">"${sample_name}"\n""^"${tag1} >> sampledata/seq2_R1_index.fasta
    echo -e ">"${sample_name}"\n""^"${tag2} >> sampledata/seq2_R2_index.fasta
  fi
  count=`expr $count + 1`
done < ${TAG_DATA2}
## ----------------------- For sequence data 3
count=0
while read row || [ -n "${row}" ]; do
  if ((count > 0)); then
    sample_name=$(echo ${row} | cut -d , -f 1)
    tag1=$(echo ${row} | cut -d , -f 2)
    tag2=$(echo ${row} | cut -d , -f 3)
    echo -e ">"${sample_name}"\n""^"${tag1} >> sampledata/seq3_R1_index.fasta
    echo -e ">"${sample_name}"\n""^"${tag2} >> sampledata/seq3_R2_index.fasta
  fi
  count=`expr $count + 1`
done < ${TAG_DATA3}
## ----------------------- For sequence data 4
count=0
while read row || [ -n "${row}" ]; do
  if ((count > 0)); then
    sample_name=$(echo ${row} | cut -d , -f 1)
    tag1=$(echo ${row} | cut -d , -f 2)
    tag2=$(echo ${row} | cut -d , -f 3)
    echo -e ">"${sample_name}"\n""^"${tag1} >> sampledata/seq4_R1_index.fasta
    echo -e ">"${sample_name}"\n""^"${tag2} >> sampledata/seq4_R2_index.fasta
  fi
  count=`expr $count + 1`
done < ${TAG_DATA4}
## ----------------------- For sequence data 5
count=0
while read row || [ -n "${row}" ]; do
  if ((count > 0)); then
    sample_name=$(echo ${row} | cut -d , -f 1)
    tag1=$(echo ${row} | cut -d , -f 2)
    tag2=$(echo ${row} | cut -d , -f 3)
    echo -e ">"${sample_name}"\n""^"${tag1} >> sampledata/seq5_R1_index.fasta
    echo -e ">"${sample_name}"\n""^"${tag2} >> sampledata/seq5_R2_index.fasta
  fi
  count=`expr $count + 1`
done < ${TAG_DATA5}


# ------------------------------------ #
# Demultiplexing
# ------------------------------------ #
# Move to working directory
cd ${WORKING_DIR}

# Do dual-unique-index demultiplexing using cutadapt
## ----------------------- For sequence data 1
cutadapt \
-j 36 \
-e 0 --no-indels \
--pair-adapters \
-g file:$(dirname ${TAG_DATA1})/seq1_R1_index.fasta \
-G file:$(dirname ${TAG_DATA1})/seq1_R2_index.fasta \
-o ${DEMILTIPLEX_DIR}/{name}_R1.fastq.gz \
-p ${DEMILTIPLEX_DIR}/{name}_R2.fastq.gz \
${SEQ_DIR}/${SEQ_FILE1_R1} ${SEQ_DIR}/${SEQ_FILE1_R2}
## ----------------------- For sequence data 2
cutadapt \
-j 36 \
-e 0 --no-indels \
--pair-adapters \
-g file:$(dirname ${TAG_DATA2})/seq2_R1_index.fasta \
-G file:$(dirname ${TAG_DATA2})/seq2_R2_index.fasta \
-o ${DEMILTIPLEX_DIR}/{name}_R1.fastq.gz \
-p ${DEMILTIPLEX_DIR}/{name}_R2.fastq.gz \
${SEQ_DIR}/${SEQ_FILE2_R1} ${SEQ_DIR}/${SEQ_FILE2_R2}
## ----------------------- For sequence data 3
cutadapt \
-j 36 \
-e 0 --no-indels \
--pair-adapters \
-g file:$(dirname ${TAG_DATA3})/seq3_R1_index.fasta \
-G file:$(dirname ${TAG_DATA3})/seq3_R2_index.fasta \
-o ${DEMILTIPLEX_DIR}/{name}_R1.fastq.gz \
-p ${DEMILTIPLEX_DIR}/{name}_R2.fastq.gz \
${SEQ_DIR}/${SEQ_FILE3_R1} ${SEQ_DIR}/${SEQ_FILE3_R2}
## ----------------------- For sequence data 4
cutadapt \
-j 36 \
-e 0 --no-indels \
--pair-adapters \
-g file:$(dirname ${TAG_DATA4})/seq4_R1_index.fasta \
-G file:$(dirname ${TAG_DATA4})/seq4_R2_index.fasta \
-o ${DEMILTIPLEX_DIR}/{name}_R1.fastq.gz \
-p ${DEMILTIPLEX_DIR}/{name}_R2.fastq.gz \
${SEQ_DIR}/${SEQ_FILE4_R1} ${SEQ_DIR}/${SEQ_FILE4_R2}
## ----------------------- For sequence data 5
cutadapt \
-j 36 \
-e 0 --no-indels \
--pair-adapters \
-g file:$(dirname ${TAG_DATA5})/seq5_R1_index.fasta \
-G file:$(dirname ${TAG_DATA5})/seq5_R2_index.fasta \
-o ${DEMILTIPLEX_DIR}/{name}_R1.fastq.gz \
-p ${DEMILTIPLEX_DIR}/{name}_R2.fastq.gz \
${SEQ_DIR}/${SEQ_FILE5_R1} ${SEQ_DIR}/${SEQ_FILE5_R2}

# Delete unknown.fastq.gz
rm seqdata_demultiplexed/unknown*.gz


# ------------------------------------------------------------------------------------- #
# Remove primers
# ------------------------------------------------------------------------------------- #
# Create trimmed data folder
mkdir seqdata_primer_trimmed
# Move to seqdata_demul
cd ${DEMILTIPLEX_DIR}

# Remove primers (allowed error rate = 10% [default])
## Forward direction
for file in *_R1.fastq.gz; do
cutadapt -j 36 \
-a ^AGTTACYYTAGGGATAACAGCG...ACRTGATCTGAGTTCAGACCGG \
-A ^CCGGTCTGAACTCAGATCAYGT...CGCTGTTATCCCTARRGTAACT \
-n 2 \
--discard-untrimmed \
-o ../seqdata_primer_trimmed/${file%_R1.fastq.gz}_R1_F.fastq.gz \
-p ../seqdata_primer_trimmed/${file%_R1.fastq.gz}_R2_F.fastq.gz \
${file} \
${file%_R1.fastq.gz}_R2.fastq.gz
done
## Reverse direction
for file in *_R1.fastq.gz; do
cutadapt -j 36 \
-a ^CCGGTCTGAACTCAGATCAYGT...CGCTGTTATCCCTARRGTAACT \
-A ^AGTTACYYTAGGGATAACAGCG...ACRTGATCTGAGTTCAGACCGG \
-n 2 \
--discard-untrimmed \
-o ../seqdata_primer_trimmed/${file%_R1.fastq.gz}_R1_R.fastq.gz \
-p ../seqdata_primer_trimmed/${file%_R1.fastq.gz}_R2_R.fastq.gz \
${file} \
${file%_R1.fastq.gz}_R2.fastq.gz
done


# ------------------------------------------------------------------------------------- #
# Merge reversed sequences
# ------------------------------------------------------------------------------------- #
# Create output directory
cd ${WORKING_DIR}
mkdir seqdata_seq_reversed
cd seqdata_primer_trimmed

# Merge reversed sequences
for file in *_R1_F.fastq.gz; do
cat ${file} ${file%R1_F.fastq.gz}R2_R.fastq.gz > ../seqdata_seq_reversed/${file%R1_F.fastq.gz}R1.fastq.gz
cat ${file%R1_F.fastq.gz}R2_F.fastq.gz ${file%R1_F.fastq.gz}R1_R.fastq.gz > ../seqdata_seq_reversed/${file%R1_F.fastq.gz}R2.fastq.gz
done


# ------------------------------------------------------------------------------------- #
# Summarize sequence reads
# ------------------------------------------------------------------------------------- #
# Move to working directory
cd ${WORKING_DIR}
mkdir 01_DemultiplexedSummary

# Summarize using seqkit
## Original sequence data
seqkit stat -a -T seqdata/*.gz > 00_DemultiplexedSummary/01_original_seq.tsv
## demultiplexed sequence data
seqkit stat -a -T seqdata_demultiplexed/*.gz > 00_DemultiplexedSummary/02_demultiplexed_seq.tsv
## primer trimmed sequence data
seqkit stat -a -T seqdata_primer_trimmed/*.gz > 00_DemultiplexedSummary/03_primertrimmed_seq.tsv
## Reversed sequences converted and merged
seqkit stat -a -T seqdata_seq_reversed/*.gz > 00_DemultiplexedSummary/04_seqreversed_seq.tsv


# ------------------------------------------------------------------------------------- #
# Delete temporal sequence files
# ------------------------------------------------------------------------------------- #
rm seqdata_demultiplexed/*.fastq.gz
rm seqdata_primer_trimmed/*.fastq.gz
