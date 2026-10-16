"""Audio EQ Cookbook LPF.

Analog prototype (Q = 1/√2 → 2nd-order Butterworth):

    H(s) = ωc² / (s² + (ωc/Q)·s + ωc²)

The cookbook discretizes it with the bilinear transform s = 2fs·(1 - z⁻¹)/(1 + z⁻¹),
with ωc prewarped so the cutoff lands exactly on f0.

Usage:
    uv run lpf                       # f0 = 1 kHz, Q = 0.707, fs = 48 kHz
    uv run lpf --f0 300 --q 4
    uv run lpf --save figs/          # save PNGs instead of opening windows
"""

import numpy as np

from .cli import Filter, run
from .rbj import FS, normalize, prewarp, rbj_params


def digital(f0: float, q: float, fs: float = FS, gain_db: float = 0.0):
    """Cookbook coefficients, normalized by a0. gain_db is unused."""
    p = rbj_params(f0, q, fs)
    c, alpha = p.cos_w0, p.alpha
    b = [(1 - c) / 2, 1 - c, (1 - c) / 2]
    a = [1 + alpha, -2 * c, 1 - alpha]
    return normalize(b, a)


def analog(f0: float, q: float, fs: float = FS, gain_db: float = 0.0):
    """Analog prototype with prewarped ωc. Returns (num, den, ωc)."""
    wc = prewarp(f0, fs)
    return np.array([wc**2]), np.array([1.0, wc / q, wc**2]), wc


FILTER = Filter(name="LPF", label="LPF · 12 dB/oct", digital=digital, analog=analog)


def main():
    run(FILTER)


if __name__ == "__main__":
    main()