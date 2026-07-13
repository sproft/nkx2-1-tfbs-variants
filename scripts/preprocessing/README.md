# EMSA-seq preprocessing

The pipeline that turns the raw sequencing reads into the log2 fold-change training labels, from the
scripts written by Freyda Lenihan-Geels and run for the study. A few minimal corrections were applied
for this release so that the scripts run as intended; they are listed at the end and none change the
analysis. The resulting DESeq2 tables are also deposited at GEO `GSE319141`, so the models can be
retrained without rerunning this pipeline.

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
matching their own design. Run the block that matches the library you are processing. Set `EMSA_ROOT`
to the directory holding the raw FASTQ files before running:

```bash
EMSA_ROOT=/path/to/raw_fastq sh EMSAseq_raw_data_analysis.sh
```

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

## Corrections applied for this release

The scripts are otherwise as deposited. The following minimal corrections were made so that they run
as intended; none change the analysis.

- **`EMSAseq_deseq2_analysis.R`, line 100:** `name0` corrected to `name=` in the `res_0.5` call,
  which was otherwise a syntax error (ALL library only).
- **`EMSAseq_deseq2_analysis.R`, lines 46 to 52:** the FLANK merge chain now accumulates into
  `dat_flank` rather than `dat`, matching the CORE and ALL blocks, so all FLANK samples are included
  rather than only the first two.
- **`EMSAseq_raw_data_analysis.sh`:** the hardcoded `/massdata/...` paths are replaced by an
  `EMSA_ROOT` variable (`rootdir="${EMSA_ROOT:?...}"`, `results="$rootdir/processed"`), and the three
  library-specific final steps now read the step-5 output from `$results` consistently.

Still done by hand in the R session, as in the original: import the per-sample count tables; assign
`dat` to the working library (`dat_core`, `dat_flank` or `dat_all`); define `coldata`; and
`library(RColorBrewer)` for the heatmap palette. The three library-specific blocks in the shell
script share one output filename, so run only the block for the library you are processing.
