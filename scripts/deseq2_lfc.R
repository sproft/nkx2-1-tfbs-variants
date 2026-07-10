# Placeholder: DESeq2 + apeglm computation of EMSA-seq log2 fold changes.
#
# This script is not yet available in this repository. It computes, per mutant library
# (CORE, FLANK, ALL), the log2 fold change in enrichment between the bound and input
# pools across the molar ratios, with apeglm shrinkage, and writes the tables that
# `preprocess.py` consumes as training labels.
#
# Expected inputs:  per-sequence counts for each library and molar ratio
# Expected outputs: <library>.csv with columns `sequence` and `log2FC`
#
# Until it lands here, the resulting tables are available from GEO accession GSE319141,
# which holds the raw reads, the sequence counts, and the DESeq2 results.
#
# Software versions used (Supplementary Table S2): DESeq2 1.50.0, apeglm 1.32.0.

stop("scripts/deseq2_lfc.R is not yet included in this repository; see GEO GSE319141 for the DESeq2 output tables.")
