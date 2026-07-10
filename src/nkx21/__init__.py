"""Trained models and code for the NKX2-1 transcription factor binding site study.

Quickstart::

    from nkx21 import load_model, predict

    flank = load_model("flank")                 # the best-generalizing model
    log2fc = predict(flank, ["CTCAAGTGCTGGTA"])  # 14 bp variable sequence
"""

from .encoding import INPUT_LENGTH, add_context, one_hot, reverse_complement
from .loader import LIBRARIES, LoadedModel, Scaling, load_baseline, load_checkpoint, load_model
from .models import ARCHITECTURES, RNN, VCNNBPNet
from .predict import predict, score_windows

__all__ = [
    "ARCHITECTURES", "INPUT_LENGTH", "LIBRARIES", "LoadedModel", "RNN", "Scaling", "VCNNBPNet",
    "add_context", "load_baseline", "load_checkpoint", "load_model", "one_hot",
    "predict", "reverse_complement", "score_windows",
]
__version__ = "1.0.0"
