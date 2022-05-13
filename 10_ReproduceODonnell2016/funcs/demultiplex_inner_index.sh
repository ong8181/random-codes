#!/bin/bash

#### Manual demultiplex script for "eDNA-early-pooling"
####

#---------------------------------------------------#
# REQUIRED: seqkit (https://bioinf.shenwei.me/seqkit)
# REQUIRED: fastp (https://github.com/OpenGene/fastp)
#---------------------------------------------------#

# For dev
#SEQ_DIR=seqdata_SingleEnd
#SEQ_FILE_R1=sample_fastq_R1_SE.fastq.gz
#SEQ_FILE_R2=sample_fastq_R2.fastq.gz
#SAMPLE_DATA=sampledata/index_info2.csv
#OUTPUT_DIR=demultiplexSE_Out
#OUTPUT_DIR=${5:-demultiplex_Out}
#SEQ_MODE=SE_SingleID
#SEQ_MODE=${6:-PE_DualID}


# Define function
function demultiplex_inner_index () {
  # Specify parameters
  #WORKING_DIR=$1
  # SEQ_MODE=SE_SingleID
  # Or, PE_SingleID, SE_SingleID, SE_DualID
  SEQ_MODE=$1
  SAMPLE_DATA=$2
  OUTPUT_DIR=$3
  SEQ_DIR=$4
  SEQ_FILE_R1=$5
  SEQ_FILE_R2=${6:-DUMMY}
  
  # Get original sequence file name
  #SEQ_FILE_R1=""
  #SEQ_FILE_R2=""
  
  # Create output directory
  mkdir ${OUTPUT_DIR}
  
  
  # ------------------------------------------------------------------------------------- #
  # Step 1. Extract barcode region of fastq (the first 1-8 bp region)
  # ------------------------------------------------------------------------------------- #
  # Prepare an output folder
  mkdir ${OUTPUT_DIR}/01_Out
  cd ${SEQ_DIR}
  
  # Manual demultiplexing only for Undetermined_*.fastq.gz
  echo -e "\nSEQ_MODE:" ${SEQ_MODE}
  echo -e "Extracting index sequences...\n"
  #for file in *.fastq.gz; do
  #  seqkit subseq -r 1:8 ${file} | gzip -c > ../${OUTPUT_DIR}/01_Out/${file%.fastq.gz}_index.fastq.gz
  #done
  if [[ "$SEQ_MODE" == "PE_DualID" ]] || [[ "$SEQ_MODE" == "SE_DualID" ]]; then
    seqkit subseq -r 1:8 ${SEQ_FILE_R1} | gzip -c > ../${OUTPUT_DIR}/01_Out/R1_index.fastq.gz
    seqkit subseq -r 1:8 ${SEQ_FILE_R2} | gzip -c > ../${OUTPUT_DIR}/01_Out/R2_index.fastq.gz
  elif  [[ "$SEQ_MODE" == "PE_SingleID" ]] || [[ "$SEQ_MODE" == "SE_SingleID" ]]; then
    seqkit subseq -r 1:8 ${SEQ_FILE_R1} | gzip -c > ../${OUTPUT_DIR}/01_Out/R1_index.fastq.gz
  fi
  
  
  # ------------------------------------------------------------------------------------- #
  # Step 2. Extract IDs of index sequences
  # ------------------------------------------------------------------------------------- #
  # Prepare an output folder
  cd ../${OUTPUT_DIR}/01_Out
  mkdir ../02_Out
  
  ### Extracting IDs
  echo "SEQ_MODE:" ${SEQ_MODE}
  echo -e "Extracting sequence IDs for each sample...\n"

  ## (One may perform quality filteirng of the ID sequences here using e.g., fastp)
  ## (Quality filtering may not be necessary for dual-unique index)
  
  #### Set counter to exclude CSV header
  count=0
  while read row || [ -n "${row}" ]; do
    if ((count > 0)); then
      sample_name=$(echo ${row} | cut -d , -f 1)
      if [[ "$SEQ_MODE" == "PE_DualID" ]] || [[ "$SEQ_MODE" == "SE_DualID" ]]; then
        index1=$(echo ${row} | cut -d , -f 2)
        index2=$(echo ${row} | cut -d , -f 3)
        seqkit grep --quiet -srip ^$index1 R1_index.fastq.gz -o ../02_Out/${sample_name}_R1_ID.fastq.gz
        seqkit grep --quiet -srip ^$index2 R2_index.fastq.gz -o ../02_Out/${sample_name}_R2_ID.fastq.gz
        #seqkit grep --quiet -srip ^$index1 R1_index.fastq.gz -o ../02_Out/Raw_${sample_name}_R1_ID.fastq.gz
        #seqkit grep --quiet -srip ^$index2 R2_index.fastq.gz -o ../02_Out/Raw_${sample_name}_R2_ID.fastq.gz
        #fastp -G -A -L -q 15 -u 37 --in1 ../02_Out/Raw_${sample_name}_R1_ID.fastq.gz --out1 ../02_Out/${sample_name}_R1_ID.fastq.gz
        #fastp -G -A -L -q 15 -u 37 --in1 ../02_Out/Raw_${sample_name}_R2_ID.fastq.gz --out1 ../02_Out/${sample_name}_R2_ID.fastq.gz
      elif  [[ "$SEQ_MODE" == "PE_SingleID" ]] || [[ "$SEQ_MODE" == "SE_SingleID" ]]; then
        index1=$(echo ${row} | cut -d , -f 2)
        seqkit grep --quiet -srip ^$index1 R1_index.fastq.gz -o ../02_Out/${sample_name}_R1_ID.fastq.gz
        #seqkit grep --quiet -srip ^$index1 R1_index.fastq.gz -o ../02_Out/Raw_${sample_name}_R1_ID.fastq.gz
        #fastp -G -A -L -q 15 -u 37 --in1 ../02_Out/Raw_${sample_name}_R1_ID.fastq.gz --out1 ../02_Out/${sample_name}_R1_ID.fastq.gz
      fi
    fi
    count=`expr $count + 1`
  done < ../../${SAMPLE_DATA}

  # Delete reports and Raw index files
  #rm ../02_Out/Raw_*
  #rm *.html
  #rm *.json

  # ------------------------------------------------------------------------------------- #
  # Step 2. Get common IDs of R1 and R2 reads
  # ------------------------------------------------------------------------------------- #
  # Prepare an output folder
  cd ../02_Out
  mkdir ../03_Out
  
  ## Extract sequence IDs
  if [[ "$SEQ_MODE" == "PE_DualID" ]] || [[ "$SEQ_MODE" == "SE_DualID" ]]; then
    echo "SEQ_MODE:" ${SEQ_MODE}
    echo -e "Extracting common sequence IDs of read 1 and read 2..."
    echo -e "(Quality filteirng of index sequences is not performed, but may be added here using e.g., fastp)"
    echo -e "(Quality filtering may not be necessary for dual-unique index)\n"
    for file in *_R1_ID.fastq.gz; do
      seqkit grep --quiet -f <(seqkit seq -ni ${file}) ${file%_R1_ID.fastq.gz}_R2_ID.fastq.gz | seqkit seq -ni > ../03_Out/${file%_R1_ID.fastq.gz}_common_ID.txt
    done
  elif  [[ "$SEQ_MODE" == "PE_SingleID" ]] || [[ "$SEQ_MODE" == "SE_SingleID" ]]; then
    echo "SEQ_MODE:" ${SEQ_MODE}
    echo -e "Extracting sequence IDs of read 1..."
    echo -e "(Quality filteirng of index sequences is not performed, but may be added here using e.g., fastp)\n"
    for file in *_R1_ID.fastq.gz; do
      seqkit seq -ni ${file} > ../03_Out/${file%_R1_ID.fastq.gz}_common_ID.txt
    done
  fi

  
  # ------------------------------------------------------------------------------------- #
  # Step 3. Get sequences based on the common IDs
  # (Trim index sequences here)
  # ------------------------------------------------------------------------------------- #
  # Prepare an output folder
  cd ../03_Out
  mkdir ../04_Out
  
  ## Extract sequences based on common IDs
  ## (Index sequences are trimmed here)
  if [[ "$SEQ_MODE" == "PE_DualID" ]] || [[ "$SEQ_MODE" == "PE_SingleID" ]]; then
    echo "SEQ_MODE:" ${SEQ_MODE}
    echo "Extracting sequences for each sample based on the common IDs..."
    echo "(Index sequences are trimmed at this step...)"
    echo -e "(This step may take time...)\n"
    for file in *_common_ID.txt; do
      seqkit grep --quiet -f ${file} ../../${SEQ_DIR}/${SEQ_FILE_R1} | seqkit subseq -r 9:-1 | gzip -c > ../04_Out/${file%_common_ID.txt}_R1.fastq.gz
      seqkit grep --quiet -f ${file} ../../${SEQ_DIR}/${SEQ_FILE_R2} | seqkit subseq -r 9:-1 | gzip -c > ../04_Out/${file%_common_ID.txt}_R2.fastq.gz
    done
  elif  [[ "$SEQ_MODE" == "SE_DualID" ]] || [[ "$SEQ_MODE" == "SE_SingleID" ]]; then
    echo "SEQ_MODE:" ${SEQ_MODE}
    echo "Extracting sequences for each sample based on the extracted IDs..."
    echo "(Index sequences are trimmed at this step...)"
    echo -e "(This step may take time...)\n"
    for file in *_common_ID.txt; do
      seqkit grep --quiet -f ${file} ../../${SEQ_DIR}/${SEQ_FILE_R1} | seqkit subseq -r 9:-1 | gzip -c > ../04_Out/${file%_common_ID.txt}_R1.fastq.gz
    done
  fi

  # ------------------------------------------------------------------------------------- #
  # Move demultiplexed fastq files and clean up folders...
  # ------------------------------------------------------------------------------------- #
  echo -e "Cleaning up temporal directories...\n"
  ## Move primer-trimmed sequences
  cd ../
  mv 04_Out/*.fastq.gz ./
  
  ## Delete temporal files
  rm -r 01_Out
  rm -r 02_Out
  rm -r 03_Out
  rm -r 04_Out
  
  cd ..
}


#---------------------------------------------------#
# Execute function
#---------------------------------------------------#
demultiplex_inner_index $1 $2 $3 $4 $5 $6
# $1=sequence_mode
# $2=index_information
# $3=Output_folder_name (optional)
# $4=sequence_data_folder
# $5=R1_fastq
# $6=R2_fasta (optional)



