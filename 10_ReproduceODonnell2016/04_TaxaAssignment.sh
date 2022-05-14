####
#### Taxa assignment
#### 2022.5.12
####

#----------- Taxa assignment using claident -----------#
#cd DADA2_OUTPUT_FOLDER

FASTA_FILE="../03_OTUClusteringOut/OTU_seqs.fa"
OUTPUT_FOLDER="04_TaxaAssignmentOut"
mkdir ${OUTPUT_FOLDER}
cd ${OUTPUT_FOLDER}

# Overall_genus
#clmakecachedb 0.9.2021.10.22
clmakecachedb --blastdb=overall_genus --numthreads=72 ${FASTA_FILE} overall_genus_cache
clidentseq --blastdb=overall_genus_cache --numthreads=72 ${FASTA_FILE} overall_genus_clidentseq
classigntax --taxdb=overall_genus --maxpopposer=0.10 --minsoratio=9 overall_genus_clidentseq overall_genus_classigntax

# Check overall_class
clmakecachedb --blastdb=overall_class --numthreads=72 ${FASTA_FILE} overall_class_cache
clidentseq --blastdb=overall_class_cache --numthreads=72 ${FASTA_FILE} overall_class_clidentseq
classigntax --taxdb=overall_class --maxpopposer=0.10 --minsoratio=9 overall_class_clidentseq overall_class_classigntax

## Merge identification results (overall_class + overall_genus)
clmergeassign --priority=descend overall_genus_classigntax overall_class_classigntax merge_classigntax

# Delete large files
rm -r overall_class_cache
rm -r overall_genus_cache
