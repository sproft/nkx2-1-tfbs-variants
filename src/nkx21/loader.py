"""Loading the released checkpoints.

The released ``.pt`` files contain only what inference needs: a ``state_dict``,
the architecture name, its hyperparameters, the label-scaling parameters, and
provenance metadata. Optimizer state has been stripped, and no reference to the
original training package is retained, so the weights load under
``torch.load(..., weights_only=True)`` with no import shim.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import torch

from .models import ARCHITECTURES

LIBRARIES = ("core", "flank", "all")


@dataclass(frozen=True)
class Scaling:
    """Label scaling applied during training (z-scoring of the log2 fold changes)."""

    method: str
    mean: float
    std: float

    def inverse(self, scores: torch.Tensor) -> torch.Tensor:
        """Map standardized model output back to log2 fold change."""
        if self.method != "standardize":
            raise ValueError(f"unsupported scaling method: {self.method}")
        return scores * self.std + self.mean


@dataclass(frozen=True)
class LoadedModel:
    model: torch.nn.Module
    scaling: Scaling
    library: str
    arch: str
    hyper_parameters: dict[str, Any]
    metadata: dict[str, Any]


def load_checkpoint(path: str | Path, device: str = "cpu") -> LoadedModel:
    """Load a released checkpoint into eval mode.

    ``weights_only=True`` is the strictest safe mode: it refuses to execute any
    pickled globals. That it succeeds is the guarantee that these weights carry
    no dependency on the original training code.
    """
    ckpt = torch.load(Path(path), map_location=device, weights_only=True)

    arch = ckpt["arch"]
    if arch not in ARCHITECTURES:
        raise ValueError(f"unknown architecture {arch!r}; expected one of {sorted(ARCHITECTURES)}")

    model = ARCHITECTURES[arch](**ckpt["hyper_parameters"])
    model.load_state_dict(ckpt["state_dict"], strict=True)
    model.eval().to(device)

    sc = ckpt["scaling"]
    return LoadedModel(
        model=model,
        scaling=Scaling(sc["scaling_method"], float(sc["mean"]), float(sc["std"])),
        library=ckpt["library"],
        arch=arch,
        hyper_parameters=dict(ckpt["hyper_parameters"]),
        metadata=dict(ckpt.get("metadata", {})),
    )


def load_model(library: str, root: str | Path = "models", device: str = "cpu") -> LoadedModel:
    """Load the VCNNBPNet model for ``library`` in ``{"core", "flank", "all"}``."""
    if library not in LIBRARIES:
        raise ValueError(f"library must be one of {LIBRARIES}, got {library!r}")
    return load_checkpoint(Path(root) / library / "vcnnbpnet.pt", device=device)


def load_baseline(library: str, root: str | Path = "models", device: str = "cpu") -> LoadedModel:
    """Load the recurrent baseline for ``library`` (fixed 24 bp input)."""
    if library not in LIBRARIES:
        raise ValueError(f"library must be one of {LIBRARIES}, got {library!r}")
    return load_checkpoint(Path(root) / "baselines" / f"rnn_{library}.pt", device=device)
