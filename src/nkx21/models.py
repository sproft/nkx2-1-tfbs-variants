"""Model definitions for the NKX2-1 transcription factor binding site study.

These are pure-PyTorch reimplementations of the two architectures used in the paper.
They reproduce the released checkpoints exactly, keyed identically to the trained
``state_dict``s, and require neither PyTorch Lightning nor any import shim.

VCNNBPNet is the model introduced in the paper. Its original implementation class
was named ``VCNNBpnet``; both spellings are accepted by :data:`ARCHITECTURES`.
"""

from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F


class VCNNBPNet(nn.Module):
    """Dilated-convolution CNN (BPNet-style) with adaptive max pooling.

    The adaptive pooling layer is what makes the input length variable: any
    sequence at least as long as the receptive field is pooled to a fixed
    ``pool_output_size`` before the dense block. The models were trained on
    24 bp sequences but can score longer windows.

    Hyperparameters were tuned per library; see ``models/<library>/config.json``.
    """

    def __init__(
        self,
        input_channels: int = 4,
        num_channels: int = 64,
        kernel_size: int = 21,
        dilations=(1, 1, 2, 4, 8),
        pool_output_size: int = 128,
        dense_sizes=(128, 64),
        classify: bool = False,
        **_ignored,
    ):
        super().__init__()
        self.classify = classify

        self.conv_layers = nn.ModuleList()
        in_ch = input_channels
        for dilation in dilations:
            padding = (kernel_size - 1) // 2 * dilation
            self.conv_layers.append(
                nn.Conv1d(in_ch, num_channels, kernel_size, padding=padding, dilation=dilation)
            )
            in_ch = num_channels

        self.pool = nn.AdaptiveMaxPool1d(pool_output_size)

        self.dense_layers = nn.ModuleList()
        fc_in = pool_output_size * num_channels
        for size in dense_sizes:
            self.dense_layers.append(nn.Linear(fc_in, size))
            fc_in = size

        self.output_layer = nn.Linear(fc_in, 1)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        for conv in self.conv_layers:
            x = F.relu(conv(x))
        x = self.pool(x)
        x = x.view(x.size(0), -1)
        for dense in self.dense_layers:
            x = F.relu(dense(x))
        out = self.output_layer(x)
        return torch.sigmoid(out) if self.classify else out


class RNN(nn.Module):
    """Convolution + GRU baseline from the architecture comparison.

    Unlike VCNNBPNet this model has a fixed input length: the dense layer is sized
    from ``input_length``, so it only accepts sequences of that length (24 bp).
    """

    def __init__(
        self,
        input_length: int = 24,
        input_channels: int = 4,
        conv_out_channels: int = 100,
        kernel_size: int = 5,
        pool_size: int = 2,
        dropout_conv: float = 0.3,
        gru_hidden_size: int = 100,
        gru_num_layers: int = 1,
        bidirectional: bool = True,
        dense_size: int = 100,
        dropout_fc: float = 0.6,
        classify: bool = False,
        **_ignored,
    ):
        super().__init__()
        self.classify = classify
        self.input_length = input_length

        self.conv1 = nn.Conv1d(input_channels, conv_out_channels, kernel_size, padding="valid")
        self.maxpool1 = nn.MaxPool1d(pool_size, stride=pool_size)
        self.dropout1 = nn.Dropout(dropout_conv)
        self.gru1 = nn.GRU(
            conv_out_channels, gru_hidden_size, num_layers=gru_num_layers,
            batch_first=True, bidirectional=bidirectional,
        )

        conv_out_len = (input_length - kernel_size + 1) // pool_size
        directions = 2 if bidirectional else 1
        self.fc1 = nn.Linear(conv_out_len * gru_hidden_size * directions, dense_size)
        self.dropout2 = nn.Dropout(dropout_fc)
        self.fc2 = nn.Linear(dense_size, 1)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = F.relu(self.conv1(x))
        x = self.maxpool1(x)
        x = self.dropout1(x)
        x = x.permute(0, 2, 1)  # (batch, seq, features) for the GRU
        x, _ = self.gru1(x)
        x = x.contiguous().view(x.size(0), -1)
        x = F.relu(self.fc1(x))
        x = self.dropout2(x)
        out = self.fc2(x)
        return torch.sigmoid(out) if self.classify else out


#: Architecture registry. Both spellings of the model name resolve to the same class.
ARCHITECTURES = {"VCNNBPNet": VCNNBPNet, "VCNNBpnet": VCNNBPNet, "RNN": RNN}
