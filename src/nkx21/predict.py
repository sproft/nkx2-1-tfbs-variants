"""Scoring sequences with the released models.

Predictions are returned in log2 fold change, the units of the EMSA-seq
measurement, by inverting the standardization applied to the training labels.
Requesting raw standardized output is possible but rarely what you want.
"""

from __future__ import annotations

from typing import Iterable, Sequence

import torch

from .encoding import INPUT_LENGTH, add_context, one_hot
from .loader import LoadedModel


@torch.no_grad()
def predict(
    loaded: LoadedModel,
    seqs: Sequence[str],
    *,
    with_context: bool = True,
    inverse_transform: bool = True,
    batch_size: int = 512,
    device: str = "cpu",
) -> torch.Tensor:
    """Score sequences and return a 1-D tensor of predictions.

    Args:
        loaded: a model from :func:`nkx21.loader.load_model`.
        seqs: DNA sequences. If ``with_context`` they are 14 bp variable
            sequences and the constant thyroglobulin context is added; otherwise
            they are passed through as-is (already the full model input).
        inverse_transform: if true (default) return log2 fold change; if false
            return the model's standardized output.
    """
    if with_context:
        seqs = [add_context(s) for s in seqs]

    out = []
    for start in range(0, len(seqs), batch_size):
        x = one_hot(seqs[start : start + batch_size]).to(device)
        out.append(loaded.model(x).squeeze(-1).cpu())
    scores = torch.cat(out) if out else torch.empty(0)

    return loaded.scaling.inverse(scores) if inverse_transform else scores


@torch.no_grad()
def score_windows(
    loaded: LoadedModel,
    sequence: str,
    *,
    window: int = INPUT_LENGTH,
    stride: int = 1,
    inverse_transform: bool = True,
    device: str = "cpu",
) -> torch.Tensor:
    """Slide a window across a long sequence and score every position.

    This is how the ChIP-seq peaks were scored: the highest-scoring window
    represents the peak, mirroring how a position weight matrix is scanned.
    """
    sequence = sequence.upper()
    if len(sequence) < window:
        raise ValueError(f"sequence shorter than the window ({len(sequence)} < {window})")
    windows = [sequence[i : i + window] for i in range(0, len(sequence) - window + 1, stride)]
    return predict(
        loaded, windows, with_context=False,
        inverse_transform=inverse_transform, device=device,
    )
