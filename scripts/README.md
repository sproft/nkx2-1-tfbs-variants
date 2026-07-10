# Analysis scripts

One script per computational result. Steps 1 to 3 regenerate the models from raw data; if you only
want to *use* the models, skip to step 4 or just `load_model()` (see the top-level README).

| Step | Script | Produces | Status |
|------|--------|----------|--------|
| 1 | *(external)* `preprocess.py` from [tfbs-mcp](https://github.com/sproft/tfbs-mcp) | training tensors | see note below |
| 2 | `deseq2_lfc.R` | log2 fold changes from sequence counts | **missing, see below** |
| 3 | *(external)* Lightning CLI from tfbs-mcp | the trained checkpoints | |
| - | `export_inference_checkpoints.py` | `models/` from the training checkpoints | ✅ |
| 4 | `figure6_chipseq_roc.py` | Figure 6 (ROC panels A-D, positional bias E-H) | ✅ |

## Reproducing the preprocessing

Steps 1 and 3 live in the companion library [tfbs-mcp](https://github.com/sproft/tfbs-mcp), pinned to
an immutable commit and installed with `pip install -e ".[pipeline]"`. The exact package versions are
recorded in [`../environment/`](../environment/README.md).

The training tensors were built with `preprocess.py` from that library:

```bash
python preprocess.py \
  --data   <library>.csv \
  --output data/tensors/<library> \
  --context --target log2FC --sequence sequence --rc
```

`--context` wraps each 14 bp variable sequence in its constant flanks (`CACTGC` … `TTGA`) to give
the 24 bp model input. The data are split 64 / 16 / 20 with `random_state=42`, and `--rc` then adds
the reverse complement of every sequence **within each split**, so no sequence shares a split with
its own reverse complement.

Hyperparameters were selected by sweeps on replicate-**median** targets and the final models trained
on replicate-**mean** targets. The train, validation and test partitions are identical between the
two, so no test sequence influenced model selection.

## Figure 6 panel order

`NEG_TFS = ["GATA1", "MYOD1", "NKX2-5", "RXRA"]` maps to panels **A, B, C, D** in that order. The
positional-bias panels are **E** = ALL, **F** = CORE, **G** = FLANK, **H** = PWM.

## Missing: `deseq2_lfc.R`

The DESeq2 and apeglm step that turns EMSA-seq sequence counts into the log2 fold changes used as
training labels is not yet in this repository. It will be added here. In the meantime the resulting
tables are available from GEO accession `GSE319141`, so the models can be retrained without it.
