####
#### Commands to demultiplex bcl to fastq
####

#----------- MiSeq data ----------#
#----------- Taxa assignment using claident -----------#
# (after picking ASV by DADA2)
#cd DADA2_OUTPUT_FOLDER

ASV_FILE="ASV_seqs.fa"

cd 03_SeqAnalysisDADA2_MiSeqOut

# Check overall_class
clmakecachedb --blastdb=overall_class --numthreads=72 $ASV_FILE ASV_overall_class_cache
clidentseq --blastdb=ASV_overall_class_cache --numthreads=72 $ASV_FILE ASV_overall_class_clidentseq
#classigntax --taxdb=overall_class ASV_overall_class_clidentseq ASV_overall_class_classigntax
classigntax --taxdb=overall_class --maxpopposer=0.05 --minsoratio=19 ASV_overall_class_clidentseq ASV_overall_class_classigntax

# Overall genus
clmakecachedb --blastdb=overall_genus --numthreads=72 $ASV_FILE ASV_overall_genus_cache
clidentseq --blastdb=ASV_overall_genus_cache --numthreads=72 $ASV_FILE ASV_overall_genus_clidentseq
#classigntax --taxdb=overall_genus ASV_overall_genus_clidentseq ASV_overall_genus_classigntax
classigntax --taxdb=overall_genus --maxpopposer=0.05 --minsoratio=19 ASV_overall_genus_clidentseq ASV_overall_genus_classigntax

# Merge identification results (overall_class + overall_genus)
clmergeassign --priority=descend ASV_overall_genus_classigntax ASV_overall_class_classigntax ASV_merge_classigntax




#----------- iSeq data ----------#
#----------- Taxa assignment using claident -----------#
# (after picking ASV by DADA2)
#cd DADA2_OUTPUT_FOLDER

ASV_FILE="ASV_seqs.fa"

cd 03_SeqAnalysisDADA2_iSeqOut

# Check overall_class
clmakecachedb --blastdb=overall_class --numthreads=72 $ASV_FILE ASV_overall_class_cache
clidentseq --blastdb=ASV_overall_class_cache --numthreads=72 $ASV_FILE ASV_overall_class_clidentseq
#classigntax --taxdb=overall_class ASV_overall_class_clidentseq ASV_overall_class_classigntax
classigntax --taxdb=overall_class --maxpopposer=0.05 --minsoratio=19 ASV_overall_class_clidentseq ASV_overall_class_classigntax

# Overall genus
clmakecachedb --blastdb=overall_genus --numthreads=72 $ASV_FILE ASV_overall_genus_cache
clidentseq --blastdb=ASV_overall_genus_cache --numthreads=72 $ASV_FILE ASV_overall_genus_clidentseq
#classigntax --taxdb=overall_genus ASV_overall_genus_clidentseq ASV_overall_genus_classigntax
classigntax --taxdb=overall_genus --maxpopposer=0.05 --minsoratio=19 ASV_overall_genus_clidentseq ASV_overall_genus_classigntax

# Merge identification results (overall_class + overall_genus)
clmergeassign --priority=descend ASV_overall_genus_classigntax ASV_overall_class_classigntax ASV_merge_classigntax

