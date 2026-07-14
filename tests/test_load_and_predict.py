"""Checks that the released weights load and score correctly.

These run on CPU in seconds and require no data beyond the repository itself.
"""

from __future__ import annotations

import json
import random
import sys
from pathlib import Path

import pytest
import torch

from nkx21 import (
    INPUT_LENGTH, add_context, load_checkpoint, load_model,
    one_hot, predict, reverse_complement, score_windows,
)

ROOT = Path(__file__).resolve().parents[1]
MODELS = ROOT / "models"
LIBRARIES = ["core", "flank", "all"]

# Range of the raw training labels per library, from the training tensors.
LABEL_RANGE = {"core": (-0.43, 3.50), "flank": (-4.82, 5.04), "all": (-4.80, 1.31)}


def random_seqs(n: int = 32, k: int = 14) -> list[str]:
    rng = random.Random(0)
    return ["".join(rng.choice("ACGT") for _ in range(k)) for _ in range(n)]


@pytest.mark.parametrize("library", LIBRARIES)
def test_weights_load_without_the_training_package(library: str) -> None:
    """weights_only=True refuses pickled globals, so this proves the weights are self-contained."""
    saved = {m: sys.modules.get(m) for m in ("classes", "tfbs")}
    sys.modules["classes"] = None  # type: ignore[assignment]
    sys.modules["tfbs"] = None  # type: ignore[assignment]
    try:
        ckpt = torch.load(MODELS / library / "vcnnbpnet.pt", map_location="cpu", weights_only=True)
    finally:
        for m, v in saved.items():
            sys.modules.pop(m, None) if v is None else sys.modules.__setitem__(m, v)
    assert ckpt["arch"] == "VCNNBpnet"
    assert ckpt["scaling"]["scaling_method"] == "standardize"


@pytest.mark.parametrize("library", LIBRARIES)
def test_architecture_matches_weights(library: str) -> None:
    """load_state_dict(strict=True) fails loudly if the vendored definition drifts."""
    loaded = load_model(library, root=MODELS)
    cfg = json.loads((MODELS / library / "config.json").read_text())
    assert loaded.hyper_parameters == cfg["hyper_parameters"]
    assert loaded.library == library


@pytest.mark.parametrize("library", LIBRARIES)
def test_predictions_lie_in_the_training_label_range(library: str) -> None:
    loaded = load_model(library, root=MODELS)
    scores = predict(loaded, random_seqs())
    assert torch.isfinite(scores).all()
    lo, hi = LABEL_RANGE[library]
    assert lo - 1.0 <= scores.min().item() and scores.max().item() <= hi + 1.0


@pytest.mark.parametrize("library", LIBRARIES)
def test_inverse_transform_is_the_documented_affine_map(library: str) -> None:
    loaded = load_model(library, root=MODELS)
    seqs = random_seqs(8)
    raw = predict(loaded, seqs, inverse_transform=False)
    log2fc = predict(loaded, seqs, inverse_transform=True)
    expected = raw * loaded.scaling.std + loaded.scaling.mean
    torch.testing.assert_close(log2fc, expected)


@pytest.mark.parametrize("library", LIBRARIES)
def test_reverse_complement_agreement(library: str) -> None:
    """Trained with reverse-complement augmentation, so a site and its RC should score alike."""
    loaded = load_model(library, root=MODELS)
    x = one_hot([add_context(s) for s in random_seqs(64)])
    with torch.no_grad():
        fwd = loaded.model(x).squeeze(-1)
        rev = loaded.model(reverse_complement(x)).squeeze(-1)
    r = torch.corrcoef(torch.stack([fwd, rev]))[0, 1].item()
    assert r > 0.75, f"reverse-complement agreement unexpectedly low for {library}: r={r:.3f}"


def test_random_sequences_do_not_bind_in_the_all_library() -> None:
    """The ALL library randomizes all 14 bp, so random sequences should score near its mean."""
    loaded = load_model("all", root=MODELS)
    mean_pred = predict(loaded, random_seqs(64)).mean().item()
    assert mean_pred < -1.5, f"expected weak binding for random sequences, got {mean_pred:.3f}"


def test_variable_input_length_is_supported() -> None:
    """VCNNBPNet pools adaptively, so it scores windows longer than the 24 bp it saw."""
    loaded = load_model("flank", root=MODELS)
    long_seq = "".join(random.Random(1).choice("ACGT") for _ in range(200))
    scores = score_windows(loaded, long_seq, window=INPUT_LENGTH)
    assert scores.shape[0] == 200 - INPUT_LENGTH + 1


def test_add_context_rejects_wrong_length_and_bad_alphabet() -> None:
    with pytest.raises(ValueError, match="14 bp variable sequence"):
        add_context("ACGT")
    with pytest.raises(ValueError, match="non-ACGT"):
        add_context("ACGTNACGTACGTA")
