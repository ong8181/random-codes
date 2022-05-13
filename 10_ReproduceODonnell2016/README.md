# Codes to reproduce and re-analyze data from O'Donnell et al (2016) _PLoS ONE_
This repository includes R and shell codes to reproduce and re-analyze data from O'Donnell et al. (2016) _PLoS ONE_ (https://doi.org/10.1371/journal.pone.0148698). Some metadata can be downloaded from https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0148698#sec009. Sequence data can be downloaded from Dryad (https://datadryad.org/stash/dataset/doi%253A10.5061%252Fdryad.mp040).

## If you want to reproduce the analysis by yourself
1. Download sequence data from https://datadryad.org/stash/dataset/doi%253A10.5061%252Fdryad.mp040.
2. Decompress the file and move the sequence data (files in `doi_10.5061_dryad.mp040__v1` folder) to `seqdata` folder in this repository.
3. Run `01_DemultiplexPrimerTrim.sh` to perform demultiplexing and primer trimming. Please change the working directory path specified in the code. `cutadapt` (https://cutadapt.readthedocs.io/) and `seqkit` (https://bioinf.shenwei.me/seqkit/) are required.
4. Run `02_DADA2_PE.R` to denoise and generate ASVs. Please change the path to the demultiplexed and primer-trimmed sequence data.
5. Run `03_OTUClustering.R` to generate 97% OTUs.
6. Run `04_TaxaAssignment.sh` to assign taxa. `Claident` (https://github.com/astanabe/Claident) is required. Running on Ubuntu is recommended.
7. Run `05_Summarize.R` and `06_ProtocolEffect.R` to visualize results.

## Notes
- The analysis procedure is not identical with that of the original one, so the results may be quantitatively different from the original ones (but I believe they are qualitatively the same).
- In general, `dada2` should be performed for each run, but all runs (maybe, five runs?) were simultaneously processed. This would probably not significantly change the results. 
- Detailed information about the packages is available in `00_SessionInfo` folder.

## Results
Figure 2, one of the main results of O'Donnell et al. (2016), was successfully reproduced. Single PCR approach introduced a large variation in "among index" treatment. You can find other results in <a href=https://github.com/ong8181/random-codes/tree/master/10_ReproduceODonnell2016/06_ProtocolEffectOut>06_ProtocolEffectOut</a>. <br>

<img src="img/Reproduce_ODonnell2016_Fig2.png" width="800px">

