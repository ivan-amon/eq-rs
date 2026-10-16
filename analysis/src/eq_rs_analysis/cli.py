"""Shared runner for every filter: CLI arguments, sanity check and plots.

Each filter module (lpf.py, hpf.py, ...) defines a `Filter` and calls `run()`.
"""

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

from . import plots
from .rbj import FS

# (f0, q, fs, gain_db) -> (b, a)
DigitalFn = Callable[[float, float, float, float], tuple[np.ndarray, np.ndarray]]
# (f0, q, fs, gain_db) -> (num, den, Ωc)
AnalogFn = Callable[[float, float, float, float], tuple[np.ndarray, np.ndarray, float]]


@dataclass(frozen=True)
class Filter:
    name: str  # "LPF": used in titles and as the PNG file prefix
    label: str  # DAW view label, e.g. "LPF · 12 dB/oct"
    digital: DigitalFn  # cookbook coefficients
    analog: AnalogFn  # prewarped analog prototype
    has_gain: bool = False  # True for peaking and shelving filters (adds --gain)


def run(flt: Filter) -> None:
    ap = argparse.ArgumentParser(description=f"Analysis of the Audio EQ Cookbook {flt.name}")
    ap.add_argument("--f0", type=float, default=1000.0, help="cutoff / center frequency (Hz)")
    ap.add_argument("--q", type=float, default=1 / np.sqrt(2), help="quality factor")
    ap.add_argument("--fs", type=float, default=FS, help="sample rate (Hz)")
    if flt.has_gain:
        ap.add_argument("--gain", type=float, default=6.0, help="gain (dB)")
    ap.add_argument("--save", type=Path, help="directory to save PNGs into instead of showing windows")
    views = ap.add_argument_group("views", "pick which plots to show (combinable; default: all)")
    views.add_argument("-a", "--analog", action="store_true", help="s-plane and analog frequency response")
    views.add_argument("-d", "--digital", action="store_true", help="z-plane and digital frequency response")
    views.add_argument("--daw", action="store_true", help="DAW-style view")
    args = ap.parse_args()
    show_all = not (args.analog or args.digital or args.daw)
    gain = args.gain if flt.has_gain else 0.0

    b, a = flt.digital(args.f0, args.q, args.fs, gain)
    num, den, wc = flt.analog(args.f0, args.q, args.fs, gain)

    # scipy's bilinear transform of the prototype must match the cookbook coefficients
    bz, az = signal.bilinear(num, den, args.fs)
    assert np.allclose(b, bz) and np.allclose(a, az), "cookbook ≠ bilinear(prototype)"

    gain_txt = f"  gain = {gain:g} dB" if flt.has_gain else ""
    print(f"{flt.name}  f0 = {args.f0:g} Hz  Q = {args.q:.4f}  fs = {args.fs:g} Hz{gain_txt}")
    print(f"  Ωc (prewarped) = {wc:.3f} rad/s  →  {wc / (2 * np.pi):.3f} Hz")
    print(f"  b = {b}")
    print(f"  a = {a}")
    print(f"  z poles = {np.roots(a)}  |p| = {np.abs(np.roots(a))}")

    # lazy: each figure is only created if its view was selected
    builders = {
        "1_s_plane": (args.analog, lambda: plots.s_plane(num, den, wc, flt.name)),
        "2_analog_fourier": (args.analog, lambda: plots.fourier_analog(num, den, wc, args.fs, flt.name)),
        "3_z_plane": (args.digital, lambda: plots.z_plane(b, a, flt.name)),
        "4_digital_freq": (args.digital, lambda: plots.freq_digital(b, a, args.fs, args.f0, flt.name)),
        "5_daw": (args.daw, lambda: plots.daw_view(b, a, args.fs, args.f0, args.q, flt.label)),
    }
    figs = {key: build() for key, (selected, build) in builders.items() if show_all or selected}

    if args.save:
        args.save.mkdir(parents=True, exist_ok=True)
        prefix = flt.name.lower()
        for key, fig in figs.items():
            fig.savefig(args.save / f"{prefix}_{key}.png", dpi=130, facecolor=fig.get_facecolor())
        print(f"Saved {len(figs)} figures to {args.save}/")
    else:
        plt.show()