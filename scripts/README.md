# Analysis scripts

One script per computational result. Steps 1 to 5 regenerate the models from raw data; if you only
want to *use* the models, skip to `load_model()` (see the top-level README).

| Step | Script | Produces |
|------|--------|----------|
| 1 | `preprocessing/EMSAseq_raw_data_analysis.sh` | processed 14 bp reads from the raw FASTQ |
| 2 | `preprocessing/EMSAseq_count.py` | per-sample unique-sequence count tables |
| 3 | `preprocessing/EMSAseq_deseq2_analysis.R` | per-ratio log2 fold changes (DESeq2 + apeglm) |
| 4 | *(external)* `preprocess.py` from [tfbs-mcp](https://github.com/sproft/tfbs-mcp) | training tensors |
| 5 | *(external)* Lightning CLI from tfbs-mcp | the trained checkpoints |
| - | `export_inference_checkpoints.py` | `models/` from the training checkpoints |
| 6 | `figure6_chipseq_roc.py` | Figure 6 (ROC panels A-D, positional bias E-H) |

Steps 1 to 3 are the EMSA-seq preprocessing pipeline; see
[`preprocessing/README.md`](preprocessing/README.md) for the sample naming, the per-library
variants, and how to run the interactive DESeq2 script. Steps 4 and 5 live in the companion library.

## Reproducing the preprocessing

Steps 1 and 3 live in [tfbs-mcp](https://github.com/sproft/tfbs-mcp), a companion library **in
preparation (not yet formally released)**, pinned to an immutable commit and installed with
`pip install -e ".[pipeline]"`. The exact package versions are recorded in
[`../environment/`](../environment/README.md).

The training tensors were built with `preprocess.py` from that library:

```bash
python preprocess.py \
  --data   <library>.csv \
  --output data/tensors/<library> \
  --context --target log2FC --sequence sequence --rc
```

`--context` wraps each 14 bp variable sequence in its constant flanks (`CACTGC` and `TTGA`) to give
the 24 bp model input. The data are split 64 / 16 / 20 with `random_state=42`, and `--rc` then adds
the reverse complement of every sequence **within each split**, so no sequence shares a split with
its own reverse complement.

Hyperparameters were selected by sweeps on replicate-**median** targets and the final models trained
on replicate-**mean** targets. The train, validation and test partitions are identical between the
two, so no test sequence influenced model selection.

## Figure 6 panel order

`NEG_TFS = ["GATA1", "MYOD1", "NKX2-5", "RXRA"]` maps to panels **A, B, C, D** in that order. The
positional-bias panels are **E** = ALL, **F** = CORE, **G** = FLANK, **H** = PWM.

## From counts to training labels

The DESeq2 step produces log2 fold changes per molar ratio. These are averaged across ratios into the
per-library `<library>.csv` files that `preprocess.py` reads (the `_mean` aggregation). The DESeq2
output tables are also deposited at GEO `GSE319141`.
