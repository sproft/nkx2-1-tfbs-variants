"""Export the training checkpoints to inference-only weights.

The Lightning checkpoints in the project workspace are roughly three times larger
than they need to be (two thirds is Adam optimizer state) and embed a reference to
the training package (``classes.nn.models.VCNNBpnet``), which forces an import shim
on anyone loading them.

This script strips both. The result loads under ``torch.load(weights_only=True)``,
the strictest safe mode, which refuses to resolve any pickled global. The exported
file additionally carries the label-scaling parameters, so predictions can be mapped
back to log2 fold change without a side-car file.

Usage::

    python scripts/export_inference_checkpoints.py --source /path/to/TFBS-NN
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import torch

LIBRARIES = {"core": "core_mean", "flank": "flank_mean", "all": "all_mean"}

# Verified on the cluster against the released checkpoints (see README "Reported metrics").
# `val` is the model-selection metric; `test` is the held-out estimate.
METRICS = {
    "VCNNBpnet": {
        "core":  {"val_mse": 0.045901, "test_mse": 0.123174, "test_pearson": 0.948, "test_spearman": 0.875},
        "flank": {"val_mse": 0.051444, "test_mse": 0.052198, "test_pearson": 0.973, "test_spearman": 0.935},
        "all":   {"val_mse": 0.284889, "test_mse": 0.291039, "test_pearson": 0.842, "test_spearman": 0.679},
    },
}

# Rows include the reverse complement of every sequence, added independently within
# each split. Unique sequences are therefore half the row count.
SPLITS = {
    "core":  {"n_train": 326,       "n_val": 82,      "n_test": 104,    "space": 4 ** 4},
    "flank": {"n_train": 1_341_838, "n_val": 335_460, "n_test": 419_326, "space": 4 ** 10},
    "all":   {"n_train": 139_154,   "n_val": 34_790,  "n_test": 43_486,  "space": 4 ** 14},
}

DROP_HPARAMS = ("_class_path", "_instantiator")


def _clean_hparams(hp: dict) -> dict:
    return {k: v for k, v in hp.items() if k not in DROP_HPARAMS}


def _coverage(lib: str) -> dict:
    s = SPLITS[lib]
    unique = (s["n_train"] + s["n_val"] + s["n_test"]) // 2
    return {
        "unique_sequences": unique,
        "sequence_space": s["space"],
        "coverage_fraction": round(unique / s["space"], 6),
        "n_train_rows": s["n_train"], "n_val_rows": s["n_val"], "n_test_rows": s["n_test"],
    }


def export(ckpt_path: Path, out_path: Path, arch: str, lib: str, scaling: dict) -> tuple[float, float]:
    ckpt = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    hp = _clean_hparams(dict(ckpt["hyper_parameters"]))

    payload = {
        "arch": arch,
        "library": lib,
        "hyper_parameters": hp,
        "state_dict": {k: v.contiguous() for k, v in ckpt["state_dict"].items()},
        "scaling": scaling,
        "metadata": {
            "study": "Quantification of the effects of single nucleotide variants in NKX2.1 transcription factor binding sites",
            "target": "log2 fold change (EMSA-seq)",
            "input_length_trained": hp.get("input_length", 24),
            "context_left": "CACTGC", "context_right": "TTGA",
            "epoch": int(ckpt.get("epoch", -1)),
            "metrics": METRICS[arch][lib],
            "data": _coverage(lib),
        },
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    torch.save(payload, out_path)

    src_mb = ckpt_path.stat().st_size / 1e6
    out_mb = out_path.stat().st_size / 1e6
    return src_mb, out_mb


def verify(out_path: Path, reference_state: dict) -> None:
    """Reload with the strictest loader, with the training packages made unimportable."""
    saved = {m: sys.modules.get(m) for m in ("classes", "tfbs")}
    sys.modules["classes"] = None  # type: ignore[assignment]
    sys.modules["tfbs"] = None  # type: ignore[assignment]
    try:
        back = torch.load(out_path, map_location="cpu", weights_only=True)
    finally:
        for m, v in saved.items():
            if v is None:
                sys.modules.pop(m, None)
            else:
                sys.modules[m] = v

    assert set(back["state_dict"]) == set(reference_state), f"key mismatch in {out_path.name}"
    for k, v in reference_state.items():
        if not torch.equal(back["state_dict"][k], v):
            raise AssertionError(f"tensor mismatch for {k} in {out_path.name}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", required=True, type=Path, help="path to the TFBS-NN workspace")
    ap.add_argument("--out", type=Path, default=Path("models"))
    args = ap.parse_args()

    sm = args.source / "saved_models_final"
    sp_dir = sm / "scaling_params"
    if not sm.is_dir():
        raise SystemExit(f"not found: {sm}")

    total_src = total_out = 0.0
    for lib, dirname in LIBRARIES.items():
        scaling = json.loads((sp_dir / f"{dirname}_standardize_VCNNBpnet.json").read_text())
        assert scaling["scaling_method"] == "standardize", lib

        for arch, out in (("VCNNBpnet", args.out / lib / "vcnnbpnet.pt"),):
            src = sm / dirname / "standardize" / f"best_{arch}.ckpt"
            ref = torch.load(src, map_location="cpu", weights_only=False)["state_dict"]
            s_mb, o_mb = export(src, out, arch, lib, scaling)
            verify(out, ref)
            total_src += s_mb
            total_out += o_mb
            print(f"  {out.as_posix():<34} {s_mb:6.2f} MB -> {o_mb:5.2f} MB   verified")

        (args.out / lib / "scaling_params.json").write_text(json.dumps(scaling, indent=2) + "\n")
        cfg = {
            "library": lib, "arch": "VCNNBPNet",
            "hyper_parameters": _clean_hparams(
                dict(torch.load(sm / dirname / "standardize" / "best_VCNNBpnet.ckpt",
                                map_location="cpu", weights_only=False)["hyper_parameters"])
            ),
            "scaling": scaling,
            "metrics": METRICS["VCNNBpnet"][lib],
            "data": _coverage(lib),
        }
        (args.out / lib / "config.json").write_text(json.dumps(cfg, indent=2) + "\n")

    print(f"\n  total {total_src:.1f} MB -> {total_out:.1f} MB  ({total_out / total_src:.0%} of original)")


if __name__ == "__main__":
    main()
