# Model card

Three convolutional models predict the change in NKX2.1 binding, in log2 fold change, caused by
variants in its binding site. One model per EMSA-seq mutant library.

## Files

```
models/
  core/   vcnnbpnet.pt  scaling_params.json  config.json
  flank/  vcnnbpnet.pt  scaling_params.json  config.json
  all/    vcnnbpnet.pt  scaling_params.json  config.json
```

Each `.pt` contains the `state_dict`, the architecture name and hyperparameters, the label-scaling
parameters, and provenance metadata. Nothing else: no optimizer state, no reference to the training
package. They load under `torch.load(..., weights_only=True)`.

## Intended use

Scoring short DNA sequences for relative NKX2.1 binding strength, and comparing a reference sequence
against its variants. The output is a **relative** measure (log2 fold change in an EMSA-seq
enrichment), not an absolute dissociation constant.

Use the **FLANK** model unless you have a specific reason not to: it generalizes best, and it is the
model used for the in vivo ChIP-seq validation.

## Architecture

Shared: five dilated 1-D convolutions (dilations 1, 1, 2, 4, 8), ReLU, adaptive max pooling, a dense
block of 128 and 64 units, and a linear output. Filters, kernel size, and pooled length were tuned
per library.

| Library | Varies | Filters | Kernel | Pooled length | Parameters |
|---------|--------|--------:|-------:|--------------:|-----------:|
| CORE  | 4 bp core | 82 | 13 | 178 | 2,231,059 |
| FLANK | 10 bp flanks | 120 | 13 | 151 | 3,083,449 |
| ALL   | all 14 bp | 82 | 42 | 167 | 2,905,099 |

Adaptive pooling is what makes the input length variable: the models were trained on 24 bp but score
longer windows.

## Inputs

A 14 bp variable sequence is placed back into its constant thyroglobulin-promoter context:

```
CACTGC | NNNNNNNNNNNNNN | TTGA          6 + 14 + 4 = 24 bp
```

One-hot encoded with channel order **A, C, G, T** as `(batch, 4, length)`. This ordering is baked
into the weights.

## Outputs and scaling

Training labels were z-scored per library, so raw model output is in standardized units. The
inverse map `log2FC = output * std + mean` is applied by default; the parameters travel inside each
checkpoint and are duplicated in `scaling_params.json`.

| Library | mean | std | raw label range |
|---------|-----:|----:|-----------------|
| CORE  |  0.10264415 | 0.58060753 | [-0.43, 3.50] |
| FLANK | -0.23446292 | 1.33945072 | [-4.82, 5.04] |
| ALL   | -2.50352097 | 1.10953736 | [-4.80, 1.31] |

## Training data and splits

Sequences were split 64 / 16 / 20 into train, validation and test with `random_state=42`, and each
split was then augmented independently with the reverse complement of every sequence. Row counts
below include that augmentation, so unique sequences are half the total.

| Library | train | val | test | unique | sequence space | coverage |
|---------|------:|----:|-----:|-------:|---------------:|---------:|
| CORE  | 326 | 82 | 104 | 256 | 4⁴ = 256 | **100%** |
| FLANK | 1,341,838 | 335,460 | 419,326 | 1,048,312 | 4¹⁰ = 1,048,576 | 99.97% |
| ALL   | 139,154 | 34,790 | 43,486 | 108,715 | 4¹⁴ = 268,435,456 | 0.04% |

## Evaluation

Hyperparameters were selected on the validation split. The test split was used at no point during
training or selection.

| Library | Test MSE | Validation MSE | Test Pearson | Test Spearman |
|---------|---------:|---------------:|-------------:|--------------:|
| CORE  | 0.1232 | 0.0459 | 0.948 | 0.875 |
| FLANK | 0.0522 | 0.0514 | 0.973 | 0.935 |
| ALL   | 0.2910 | 0.2849 | 0.842 | 0.679 |

VCNNBPNet was selected by the lowest mean validation loss across the three libraries, from a
comparison of several architecture families.

## Limitations

- The CORE model is trained on 163 unique sequences and evaluated on 52. Its test error is estimated
  from few sequences and its validation-to-test gap is large. Prefer FLANK.
- Predictions are relative, in the units of an EMSA-seq enrichment, and are calibrated to the
  competitive binding regime of that assay. They correlate with competitive measurements such as
  selective chromatography but not with non-competitive microscale thermophoresis.
- The models are specific to NKX2.1 and do not transfer to other transcription factors.
- Scores are informative about binding, not about pathogenicity. An evaluation on non-coding ClinVar
  variants overlapping NKX2.1 peaks showed no discrimination, because those variants are mostly
  pathogenic through mechanisms other than NKX2.1 binding.
