# Data

Nothing large is committed to this repository. Everything below is either deposited in a public
archive or downloadable from its source.

## Deposited for this study

| Data | Accession |
|------|-----------|
| EMSA-seq raw reads, sequence counts, DESeq2 results (CORE, FLANK, ALL) | GEO `GSE319141` |
| Crystal structure, wildtype NKX2-1 homeodomain with DNA | PDB `9U18` |
| Crystal structure, CACG variant | PDB `9U19` |

The DESeq2 output tables in `GSE319141` are the training labels. Together with `preprocess.py` from
[tfbs-mcp](https://github.com/sproft/tfbs-mcp) they are sufficient to rebuild the training tensors
and retrain the models from scratch.

## Third-party resources

| Resource | Where | Used for |
|----------|-------|----------|
| ChIP-seq peaks, `Oth.ALL.50.<TF>.AllCell_no_overlaps.bed` | [ChIP-Atlas](https://chip-atlas.org/) | Figure 6 |
| NKX2-1 position weight matrix `MA1994.1` | [JASPAR](https://jaspar.elixir.no/matrix/MA1994.1/) | PWM baseline |
| Human reference genome, hg38 | UCSC or Ensembl | extracting peak sequences |

ChIP-Atlas peaks are the high-confidence set: MACS2 *Q* < 1e-50, pooled across cell types
(`AllCell`), with overlapping peaks removed within each file. For Figure 6 the positive set is
NKX2-1 and the negatives are GATA1, MYOD1, NKX2-5 and RXRA. Peaks are centred on their midpoint and
a symmetric 500 bp window is extracted. Regions overlapping between the positive and negative sets
are removed before scoring.

`MA1994.1` is the mouse Nkx2.1 matrix. The NKX2-1 DNA-binding domain is identical between mouse and
human, so it is a fair baseline for the human peaks.

## Layout

Downloaded files are gitignored. The scripts expect:

```
data/
  motifs/MA1994.1.meme          # small, may be committed
  chipseq/Oth.ALL.50.<TF>.AllCell_no_overlaps.bed
  genome/hg38.fa
```
