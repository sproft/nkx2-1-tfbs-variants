"""Sequence encoding, exactly as used to build the training tensors.

The EMSA-seq libraries vary a 14 bp window. Before encoding, each variable
sequence is placed back into its constant thyroglobulin-promoter context,
giving the 24 bp input the models were trained on:

    CACTGC | NNNNNNNNNNNNNN | TTGA
     6 bp        14 bp        4 bp
"""

from __future__ import annotations

from typing import Iterable, Sequence

import torch
import torch.nn.functional as F

LEFT_CONTEXT = "CACTGC"
RIGHT_CONTEXT = "TTGA"
VARIABLE_LENGTH = 14
INPUT_LENGTH = len(LEFT_CONTEXT) + VARIABLE_LENGTH + len(RIGHT_CONTEXT)  # 24

#: Channel order of the one-hot encoding. Do not change: it is baked into the weights.
NUC_TO_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}


def add_context(seq: str) -> str:
    """Wrap a 14 bp variable sequence in its constant flanking context."""
    seq = seq.upper()
    out = LEFT_CONTEXT + seq + RIGHT_CONTEXT
    invalid = set(out) - set("ACGT")
    if invalid:
        raise ValueError(f"sequence contains non-ACGT characters: {sorted(invalid)}")
    if len(out) != INPUT_LENGTH:
        raise ValueError(f"expected a {VARIABLE_LENGTH} bp variable sequence, got {len(seq)} bp")
    return out


def one_hot(seqs: Iterable[str]) -> torch.Tensor:
    """One-hot encode DNA sequences to a ``(batch, 4, length)`` float tensor."""
    seqs = [s.upper() for s in seqs]
    if len({len(s) for s in seqs}) != 1:
        raise ValueError("all sequences must have the same length")
    try:
        idx = torch.tensor([[NUC_TO_IDX[c] for c in s] for s in seqs], dtype=torch.long)
    except KeyError as exc:  # pragma: no cover - defensive
        raise ValueError(f"non-ACGT character in sequence: {exc}") from exc
    return F.one_hot(idx, num_classes=4).float().permute(0, 2, 1)


def reverse_complement(x: torch.Tensor) -> torch.Tensor:
    """Reverse-complement a one-hot batch of shape ``(batch, 4, length)``.

    Reversing the A,C,G,T channel axis maps A->T and C->G, so flipping both the
    channel and the length axis yields the reverse complement. This matches the
    augmentation applied during preprocessing, which was performed independently
    within each of the train, validation and test splits.
    """
    return torch.flip(x, [1, 2])


def decode(x: torch.Tensor) -> list[str]:
    """Inverse of :func:`one_hot` (argmax over channels). Useful for tests."""
    inv = "ACGT"
    return ["".join(inv[i] for i in row.argmax(0).tolist()) for row in x]
