# rDoAMP v0.1.0
Functions to extract amplicons from target sequences using a user-specified primer set.

(Convenient wrapper functions for `seqkit amplicon`)

# `doamp_auto()`
- Download sequences using `rentrez` package of R.
- Allow random sampling from searched sequence IDs.
- Extract amplicons using `seqkit amplicon`.
- Automatically expand degenerate primers and create primer list to allow `--max-mismatch` option of `seqkit amplicon`.
