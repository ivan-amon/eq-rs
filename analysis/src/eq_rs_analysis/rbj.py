"""Shared helpers for the Audio EQ Cookbook (Robert Bristow-Johnson)."""

from dataclasses import dataclass

import numpy as np

FS = 48_000.0


@dataclass(frozen=True)
class RbjParams:
    w0: float  # Ωo digital center / cutoff frequency (rad/sample)
    cos_w0: float
    sin_w0: float
    alpha: float  # sin(w0) / (2Q)
    A: float  # 10^(gain_db/40), only used by peaking and shelving filters


def rbj_params(f0: float, q: float, fs: float = FS, gain_db: float = 0.0) -> RbjParams:
    w0 = 2 * np.pi * f0 / fs # Ωo (rad)
    return RbjParams(
        w0=w0,
        cos_w0=np.cos(w0),
        sin_w0=np.sin(w0),
        alpha=np.sin(w0) / (2 * q),
        A=10 ** (gain_db / 40),
    )


def prewarp(f0: float, fs: float = FS) -> float:
    """Analog Ωc (rad/s) that the bilinear transform s = 2fs·(1 - z⁻¹)/(1 + z⁻¹) maps exactly onto f0."""
    return 2 * fs * np.tan(np.pi * f0 / fs)


def normalize(b, a):
    """Divide by a0 so that a[0] = 1."""
    b, a = np.asarray(b, dtype=float), np.asarray(a, dtype=float)
    return b / a[0], a / a[0]