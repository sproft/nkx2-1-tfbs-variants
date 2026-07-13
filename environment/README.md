# Exact environment

A snapshot of the environment in which the models were trained and the analyses were run, captured
from the BIH cluster on 2026-07-10.

**You do not need this to use the models.** Loading a checkpoint and scoring sequences needs only
PyTorch: `pip install -e .` from the repository root is enough. This directory exists so that the
paper's analyses can be reproduced against the exact package versions that produced them.

## Platform

The pins describe **Linux x86_64 with CUDA 12.1**, Python 3.11.10, PyTorch 2.4.0+cu121, built from a
conda environment. `requirements.txt` therefore includes twelve `nvidia-*` packages and `triton`,
which publish **no wheels for macOS or Windows**. Do not install it on those platforms; use the
plain `pip install -e .` route instead.

## Files

| File | What it is |
|------|------------|
| `requirements.txt` | 253 pinned packages, installable with `pip install -r`. Every conda-built package was resolved to a real version. |
| `requirements-pip-freeze.txt` | The raw `pip freeze`, kept for audit. 121 of its lines are conda `@ file:///...` paths and it is **not** installable directly. |
| `environment.yml` | Full conda export, including channels and conda-level packages. |
| `SHA256SUMS` | Checksums for the three files above. |

The raw `pip freeze` (`requirements-pip-freeze.txt`) is shipped exactly as the cluster produced it.
The installable `requirements.txt` is the same file with a single em-dash removed from its comment
header for house style; `SHA256SUMS` records the checksums of both as shipped here. They were named
`requirements_clean.txt` and `requirements.txt` on the cluster, and their headers still name the
`TFBS-NN` workspace they were generated from.

Verify before use:

```bash
cd environment && sha256sum -c SHA256SUMS
```

## Install (Linux, CUDA 12.1)

```bash
python -m venv .venv && source .venv/bin/activate

# 1. PyTorch first. The +cu121 local version tag matters: plain PyPI does not carry that build.
pip install torch==2.4.0+cu121 --index-url https://download.pytorch.org/whl/cu121

# 2. Everything else.
pip install -r environment/requirements.txt

# 3. This package.
pip install -e .
```

`requirements.txt` also pins `torch==2.4.0`. That is deliberate and harmless: under PEP 440 a
version specifier without a local label matches any local version, so the `+cu121` build installed
in step 1 satisfies it and pip will not replace it.

## A note on `tfbs-mcp`

The freeze contains an editable install of the `TFBS-NN` workspace at commit `3747683` and **no
`tfbs-mcp`**. That commit is the one that split the library out of the workspace, and the
environment was never reinstalled afterwards. The analyses therefore ran against the model code that
still lived inside `TFBS-NN`, under the `classes` package, rather than against `tfbs-mcp`.

The library now lives at [sproft/tfbs-mcp](https://github.com/sproft/tfbs-mcp), a companion package
that is in preparation and not yet formally released, and provides the same code under the `tfbs`
package. To rerun the preprocessing, training or attribution steps today, install it at the pinned
commit:

```bash
pip install -e ".[pipeline]"
```

This repository does not depend on it for inference. The model definitions needed to load and score
the released weights are vendored in `src/nkx21/models.py`, which is why the weights load with
`torch.load(weights_only=True)` and no import shim.
