# rDoAMP v0.1.0
R functions to extract amplicons from target sequences using a user-specified primer set.

(Convenient wrapper functions for `seqkit amplicon`)

## Prerequisite
- `seqkit` (https://bioinf.shenwei.me/seqkit/)
- `rentrez` (https://github.com/ropensci/rentrez)

## `doamp_auto()`
- Download sequences using `rentrez` package of R.
- Allow random sampling from searched sequence IDs.
- Extract amplicons using `seqkit amplicon`.
- Automatically expand degenerate primers and create a list of primer combinations to use `--max-mismatch` option of `seqkit amplicon` for degenerated primers.

#### Arguments
```
doamp_auto(search_query,
           F_primer,
           R_primer,
           n_retmax = 20,
           n_mismatch = 0,
           output_dir = "rDoAMP_Out",
           random_sampling = TRUE,
           random_sampling_seed = 1234,
           n_retidmax = n_retmax * 10,
           save_parameter = TRUE,
           save_stat = TRUE,
           overwrite_output_dir = FALSE)
```

#### Basic usage
```
doamp_auto("Trachurus AND mitochondrion AND 1000:20000[SLEN] ",
           F_primer = "GTCGGTAAAACTCGTGCCAGC",                    # MiFish-U-F
           R_primer = "CATAGTGGGGTATCTAATCCCAGTTTG",              # MiFish-U-R
           n_mismatch = 3)
```

## `doamp_custom()`
- Extract amplicons from a user-specified, custom FASTA file using `seqkit amplicon`.
- Automatically expand degenerate primers and create a list of primer combinations to use `--max-mismatch` option of `seqkit amplicon` for degenerated primers.

#### Arguments
```
doamp_custom(target_fasta,
             F_primer,
             R_primer,
             n_mismatch = 0,
             output_dir = "rDoAMP_Out",
             save_parameter = TRUE,
             save_stat = TRUE,
             overwrite_output_dir = FALSE) 
```

#### Basic usage
```
doamp_custom("YOUR_FASTA.fasta",
             F_primer = "GTCGGTAAAACTCGTGCCAGC",
             R_primer = "CATAGTGGGGTATCTAATCCCAGTTTG",
             n_mismatch = 3)
```

## `popular_primer_set`
- A list of popular primer sets.
