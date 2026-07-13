# EMSA-seq preprocessing

The pipeline that turns the raw sequencing reads into the log2 fold-change training labels. These are
the scripts as written by Freyda Lenihan-Geels and run for the study; they are provided as the record
of what produced the data. The resulting DESeq2 tables are also deposited at GEO `GSE319141`, so the
models can be retrained without rerunning this pipeline.

## Pipeline

```
raw *.fastq.gz
   |  EMSAseq_raw_data_analysis.sh   seqtk / fastq-filter / cutadapt
   v
processed 14 bp reads, one file per sample
   |  EMSAseq_count.py               count unique 14 bp sequences
   v
<sample>.txt   (sequence <tab> count)
   |  EMSAseq_deseq2_analysis.R      DESeq2 + apeglm, per molar ratio
   v
per-ratio log2 fold changes  ->  averaged across ratios  ->  <library>.csv
   |  preprocess.py (tfbs-mcp)       one-hot, context, split, reverse complement
   v
training tensors
```

The averaging of the per-ratio log2 fold changes into the `<library>.csv` files that
`preprocess.py` consumes (the `_mean` aggregation) is a downstream merge, not part of these three
scripts.

## Sample naming

Each sample is `<library><ratio><replicate>`:

| Field | Values |
|-------|--------|
| library | `A` = ALL, `C` = CORE, `F` = FLANK |
| ratio | `D` = unbound input pool; `1` = 0.1, `2` = 0.2, `5` = 0.5 protein-to-DNA molar ratio (0.5 is ALL only) |
| replicate | `1`, `2`, `3` |

So `FD2` is the FLANK unbound pool, replicate 2; `A51` is the ALL library at a 0.5 ratio, replicate 1.

## The three scripts

**`EMSAseq_raw_data_analysis.sh`** trims reads to 40 bp, drops reads with mean Phred below 30,
requires the correct reverse primer with cutadapt (removing off-target primer binding), trims to the
14 bp variable region, and keeps only reads that are exactly 14 bp. The final step is
library-specific and there are three variants in the file, one per library: ALL removes sequences
contaminating from the CORE design (`CCAGT...TGTTC`), while FLANK and CORE retain only sequences
matching their own design. Run the block that matches the library you are processing.

**`EMSAseq_count.py`** counts unique 14 bp sequences in one processed FASTQ and writes a
`sequence <tab> count` table. Run once per sample:

```bash
python EMSAseq_count.py sample.14bp.fastq.gz -o sample.txt
```

**`EMSAseq_deseq2_analysis.R`** runs DESeq2 with apeglm shrinkage, per library, comparing each bound
molar ratio against the unbound pool, and writes the shrunken log2 fold changes. It also produces the
MA plots, the PCA, and the sample-to-sample distance heatmap used in the supplement.

## Running the R script

It is an interactive analysis script, not a turnkey program. Before the DESeq2 call you must, in the
R session:

1. Import each `<sample>.txt` count table into a data frame named after the sample (`C11`, `CD1`,
   `F21`, `A51`, and so on).
2. After the per-library merge, assign the working frame: `dat <- dat_core` (or `dat_flank`,
   `dat_all`).
3. Define `coldata`, a data frame whose rows match the columns of `dat` and that carries the
   `Condition` and `Replicate` factors.
4. `library(RColorBrewer)` for the heatmap palette, in addition to the libraries the script loads.

Use `varianceStabilizingTransformation` for the CORE library and `rlog` for ALL and FLANK, as the
script notes.

## Software

seqtk, fastq-filter and cutadapt for read processing; DESeq2 1.50.0 and apeglm 1.32.0 for the
differential analysis, with pheatmap and RColorBrewer for the plots. Exact versions are in
[`../../environment/`](../../environment/README.md).

## Notes on the deposited scripts

These are provided as run. A reproducer should be aware of the following, all of which are easy to
work around in the interactive session:

- **`EMSAseq_deseq2_analysis.R`, line 100** reads `results(dds, name0"Condition_0.5_vs_Unbound")`.
  `name0` is a typo for `name=`; correct it before running. This line is for the ALL library only.
- **`EMSAseq_deseq2_analysis.R`, lines 46 to 52** (the FLANK merge chain) accumulate into `dat`
  rather than `dat_flank`, so only `FD1` and `FD2` reach `dat_flank`. Chain these merges into
  `dat_flank`, as the CORE and ALL blocks correctly do, to include the remaining FLANK samples.
- The script relies on objects set up by hand in the session (`coldata`, the per-sample count
  frames, the `dat` assignment) and on `library(RColorBrewer)`; these are not created by the script
  itself.
- **`EMSAseq_raw_data_analysis.sh`** uses absolute `/massdata/...` paths and its three
  library-specific final blocks share one output filename; set `rootdir`/`results` for your system
  and run only the block for the library you are processing.
