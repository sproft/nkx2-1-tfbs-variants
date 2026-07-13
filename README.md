# NKX2-1 transcription factor binding site variants

Trained models and code for **"Quantification of the effects of single nucleotide variants in NKX2.1 transcription factor binding sites."**

VCNNBPNet is a convolutional network that predicts how a single nucleotide variant changes the binding affinity of the transcription factor NKX2-1, trained on EMSA-seq libraries that measure binding to hundreds of thousands of variant sequences at once. This repository holds the three trained models, the code to run them, and the scripts that reproduce the computational figures.

## Install

```bash
git clone https://github.com/sproft/nkx2-1-tfbs-variants
cd nkx2-1-tfbs-variants
pip install -e .                   # inference needs only PyTorch
pip install -e ".[figures,dev]"    # to reproduce figures and run the tests
pip install -e ".[pipeline]"       # to rerun preprocessing, training or attribution
```

The `pipeline` extra pulls in [tfbs-mcp](https://github.com/sproft/tfbs-mcp), a companion library that is **in preparation and not yet formally released**, at a pinned commit. It is optional and is not needed to load or score the models; a released, versioned form will supersede the pin once it is available. To recreate the exact package versions the paper was produced with, see [`environment/`](environment/README.md).

## Quickstart

The models take a 14 bp variable sequence, place it back into its constant thyroglobulin-promoter context, and return the predicted **log2 fold change** in binding.

```python
from nkx21 import load_model, predict

flank = load_model("flank")   # the best-generalizing of the three models

predict(flank, ["CTCAAGTGCTGGTA",    # intact CAAG core
                "CTCACGTGCTGGTA"])   # CAAG -> CACG
# tensor([ 3.181, -0.333])
```

The model scores the disrupted core 3.5 log2 fold changes lower. That variant is the one whose crystal structure (PDB `9U19`) shows a rearranged protein-DNA interface, so the prediction recovers a structurally validated loss of binding from sequence alone.

`VCNNBPNet` pools adaptively and therefore accepts sequences longer than the 24 bp it was trained on, which is how ChIP-seq peaks are scored:

```python
from nkx21 import score_windows
scores = score_windows(flank, peak_sequence)   # slide a 24 bp window, score every position
```

## The models

Each library varies a different part of the binding site, and each model has its own tuned hyperparameters.

| Library | What varies | Filters | Kernel | Pooled length | Parameters |
|---------|-------------|--------:|-------:|--------------:|-----------:|
| `core`  | the 4 bp core | 82 | 13 | 178 | 2,231,059 |
| `flank` | the 10 bp flanks | 120 | 13 | 151 | 3,083,449 |
| `all`   | all 14 bp | 82 | 42 | 167 | 2,905,099 |

All three share five dilated convolutions (dilations 1, 1, 2, 4, 8), a dense block of 128 and 64 units, and a 24 bp training input. Recurrent baselines from the architecture comparison are in `models/baselines/`. See [`models/README.md`](models/README.md) for the full model card.

## Reported metrics

Hyperparameters were selected on the **validation** split; the numbers below on the **held-out test** split, which was used at no point during training or selection.

| Library | Test MSE | Validation MSE | Test Pearson | Test Spearman | Sequence-space coverage |
|---------|---------:|---------------:|-------------:|--------------:|------------------------:|
| `core`  | 0.1232 | 0.0459 | 0.948 | 0.875 | 256 of 4⁴ (**100%**) |
| `flank` | 0.0522 | 0.0514 | 0.973 | 0.935 | 1,048,312 of 4¹⁰ (99.97%) |
| `all`   | 0.2910 | 0.2849 | 0.842 | 0.679 | 108,715 of 4¹⁴ (0.04%) |

The FLANK model generalizes best. It saw essentially its entire sequence space, which is why 10 bp of flank diversity is sufficient to predict relative NKX2-1 binding strength and the far larger ALL library is not worth its sequencing cost. The CORE model's validation-to-test gap follows from its training split of only 163 unique sequences.

## Data

Nothing large is shipped here. See [`data/README.md`](data/README.md) for provenance and download instructions.

| What | Where |
|------|-------|
| EMSA-seq reads, counts, DESeq2 results | GEO [`GSE319141`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE319141) |
| Crystal structures (wildtype, CACG variant) | PDB [`9U18`](https://www.rcsb.org/structure/9U18), [`9U19`](https://www.rcsb.org/structure/9U19) |
| ChIP-seq peaks (NKX2-1 + negative controls) | ChIP-Atlas, `Oth.ALL.50.<TF>.AllCell` (MACS2 *Q* < 1e-50) |
| NKX2-1 position weight matrix | JASPAR `MA1994.1` |

## Reproducing the analyses

`scripts/` maps one script to one figure; see [`scripts/README.md`](scripts/README.md).

## How the weights are packaged

The training checkpoints were PyTorch Lightning files: two thirds of each was Adam optimizer state, and each embedded a reference to the training package, which forced an import shim on anyone loading them. The released weights are stripped to a `state_dict`, the architecture and its hyperparameters, and the label-scaling parameters. They load under `torch.load(..., weights_only=True)`, the strictest safe mode, which refuses to resolve any pickled global:

```bash
python scripts/export_inference_checkpoints.py --source /path/to/TFBS-NN   # 108 MB -> 36 MB
pytest                                                                     # verifies load + scoring
```

Predictions are returned in log2 fold change because the scaling parameters travel with the weights. Requesting `inverse_transform=False` returns the standardized output instead.

## Citation

If you use these models, please cite the paper and this archive. See [`CITATION.cff`](CITATION.cff).

## License

Everything here, code and trained model weights alike, is released under the MIT license ([`LICENSE`](LICENSE)). The weights were trained on EMSA-seq data generated for this study and deposited at GEO `GSE319141`, so no upstream license restricts them.

One third-party file is redistributed: `data/motifs/MA1994.1.meme`, from JASPAR under CC BY 4.0. Its attribution requirement is independent of this repository's license; see [`data/README.md`](data/README.md). ChIP-Atlas peaks and the reference genome are downloaded rather than redistributed, and remain under the terms of their sources.
