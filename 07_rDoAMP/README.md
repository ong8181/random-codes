# rDoAMP v0.1.0
Functions to extract amplicons from target sequences using a user-specified primer set.

(Convenient wrapper functions for `seqkit amplicon`)

# `doamp_auto()`
- Download sequences using `rentrez` package of R.
- Allow random sampling from searched sequence IDs.
- Extract amplicons using `seqkit amplicon`.
- Automatically expand degenerate primers and create a list of primer combinations to allow `--max-mismatch` option of `seqkit amplicon`.

# `doamp_custom()`
- Extract amplicons from a user-specified, custom FASTA file using `seqkit amplicon`.
- Automatically expand degenerate primers and create a list of primer combinations to allow `--max-mismatch` option of `seqkit amplicon`.

# `popular_primer_set`
- A list of popular primer sets.
