####
#### Commands to demultiplex bcl to fastq
#### The number of OTUs is ca. 5500, so skip overall_class analysis
####

#----------- Taxa assignment using claident -----------#
# (after picking ASV by DADA2)
#cd DADA2_OUTPUT_FOLDER
cd ~/Desktop/20220513_ODonnell

FASTA_FILE="../03_OTUClusteringOut/OTU_seqs.fa"
OUTPUT_FOLDER="04_TaxaAssignmentOut"
mkdir ${OUTPUT_FOLDER}
cd ${OUTPUT_FOLDER}

# Overall_genus
#clmakecachedb 0.9.2021.10.22
#clmakecachedb blastn -max_target_seqs 5000 end --blastdb=overall_genus --numthreads=36 ${FASTA_FILE} overall_genus_cache
clmakecachedb --blastdb=overall_genus --numthreads=72 ${FASTA_FILE} overall_genus_cache
clidentseq --blastdb=overall_genus_cache --numthreads=72 ${FASTA_FILE} overall_genus_clidentseq
classigntax --taxdb=overall_genus --maxpopposer=0.10 --minsoratio=9 overall_genus_clidentseq overall_genus_classigntax

# Check overall_class
# If the number of OTUs is < 3000, perform below
clmakecachedb --blastdb=overall_class --numthreads=72 ${FASTA_FILE} overall_class_cache
clidentseq --blastdb=overall_class_cache --numthreads=72 ${FASTA_FILE} overall_class_clidentseq
classigntax --taxdb=overall_class --maxpopposer=0.10 --minsoratio=9 overall_class_clidentseq overall_class_classigntax
#----------- Alternative codes ------------#
##clmakecachedb blastn -max_target_seqs 10000 end --blastdb=overall_class --numthreads=36 ${FASTA_FILE} overall_class_cache
##classigntax --taxdb=overall_class overall_class_clidentseq overall_class_classigntax
##classigntax --taxdb=overall_class --maxpopposer=0.05 --minsoratio=19 overall_class_clidentseq overall_class_classigntax
#------------------------------------------#

## Merge identification results (overall_class + overall_genus)
clmergeassign --priority=descend overall_genus_classigntax overall_class_classigntax merge_classigntax

# Delete large files
rm -r overall_class_cache
rm -r overall_genus_cache

# Move file
#mv overall_class_clidentseq ../${OUTPUT_FOLDER}
#mv overall_class_classigntax ../${OUTPUT_FOLDER}
#mv overall_genus_clidentseq ../${OUTPUT_FOLDER}
#mv overall_genus_classigntax ../${OUTPUT_FOLDER}
#mv merge_classigntax ../${OUTPUT_FOLDER}


