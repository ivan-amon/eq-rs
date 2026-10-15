"""Generic plots for an analog system H(s) and its digital counterpart H(z).

Conventions:
- Analog:  num, den are polynomials in s (descending powers), as in scipy.signal.freqs.
- Digital: b, a are coefficients in z⁻¹ (b0 + b1 z⁻¹ + ...), as in scipy.signal.freqz.
- The system is assumed causal, so the ROC is Re{s} > max Re{p} in s and |z| > max |p| in z.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
from scipy import signal

DB_FLOOR = -60.0
DB_CEIL = 30.0


# ------------------------------------------------------------------- helpers

def _db(x):
    return 20 * np.log10(np.maximum(np.abs(x), 1e-12))


def _clip(db):
    return np.clip(db, DB_FLOOR, DB_CEIL)


def _H_s(num, den, s):
    return np.polyval(num, s) / np.polyval(den, s)


def _hz_label(f, _pos=None):
    return f"{f / 1000:g}k" if abs(f) >= 1000 else f"{f:g}"


def _plot_roots(ax, roots, scale, marker, color, label):
    """Plot roots and annotate multiplicity (e.g. the double zero at z = -1)."""
    if len(roots) == 0:
        return
    pts = np.round(roots / scale, 6)
    uniq, counts = np.unique(pts, return_counts=True)
    style = dict(s=90, linewidths=2, zorder=5, label=label)
    if marker == "o":
        ax.scatter(uniq.real, uniq.imag, marker="o", facecolors="none", edgecolors=color, **style)
    else:
        ax.scatter(uniq.real, uniq.imag, marker="x", c=color, **style)
    for r, n in zip(uniq, counts):
        if n > 1:
            ax.annotate(f"×{n}", (r.real, r.imag), xytext=(8, 8), textcoords="offset points", color=color)


# ------------------------------------------------- 1. s-plane: poles, zeros, ROC

def s_plane(num, den, wc, title):
    z, p, _ = signal.tf2zpk(num, den)
    lim = 2.0

    fig, ax = plt.subplots(figsize=(7, 7))

    # |H(s)| contour lines as background
    g = np.linspace(-lim, lim, 400)
    X, Y = np.meshgrid(g, g)
    mag = _clip(_db(_H_s(num, den, (X + 1j * Y) * wc)))
    cs = ax.contour(X, Y, mag, levels=np.arange(-40, 31, 6), cmap="Greys", linewidths=0.6, alpha=0.6)
    ax.clabel(cs, fontsize=7, fmt="%d dB")

    # ROC (causal system): Re{s} > max Re{p}
    sigma_roc = p.real.max() / wc
    ax.axvspan(sigma_roc, lim, color="tab:green", alpha=0.12,
               label=f"ROC: Re{{s}} > {sigma_roc:.3f}·Ωc")
    ax.axvline(sigma_roc, color="tab:green", ls=":", lw=1.2)

    ax.axhline(0, color="k", lw=0.8)
    ax.axvline(0, color="k", lw=1.6, label="jΩ axis (Fourier)")
    ax.add_patch(Circle((0, 0), 1, fill=False, ls="--", color="gray", lw=0.8))  # |s| = Ωc

    _plot_roots(ax, p, wc, "x", "crimson", "poles")
    _plot_roots(ax, z, wc, "o", "royalblue", "zeros")

    n_inf = len(p) - len(z)
    if n_inf > 0:
        ax.text(0.02, 0.02, f"{n_inf} zero(s) at s → ∞", transform=ax.transAxes, color="royalblue")

    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_xlabel("σ / Ωc")
    ax.set_ylabel("Ω / Ωc")
    ax.set_title(f"{title} — s-plane")
    ax.legend(loc="upper right", fontsize=8)
    return fig


# ------------------------------------------------- 2. Analog Fourier H(jΩ)

def fourier_analog(num, den, wc, fs, title):
    # Same span as the digital plot: Ω ∈ [−π·fs, π·fs] ↔ ω ∈ [−π, π] (ω = Ω·Ts)
    W = np.pi * fs
    om = np.linspace(-W, W, 8193)
    _, h = signal.freqs(num, den, worN=om)
    mag = np.abs(h)
    _, hc = signal.freqs(num, den, worN=[wc])
    mc = np.abs(hc[0])

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(om, mag)
    ax.axvline(wc, color="crimson", ls="--", lw=1, label=f"±Ωc = ±{wc:.1f} rad/s")
    ax.axvline(-wc, color="crimson", ls="--", lw=1)
    ax.scatter([-wc, wc], [mc, mc], color="crimson", zorder=5)
    ax.annotate(f"|H(jΩc)| = {mc:.3f}", (wc, mc), xytext=(10, 6), textcoords="offset points",
                color="crimson")

    ax.set_xlim(-W, W)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _p: f"{x / 1000:g}k" if x else "0"))
    ax.grid(True, alpha=0.3)

    sec = ax.secondary_xaxis("top", functions=(lambda x: x / (2 * np.pi), lambda f: f * 2 * np.pi))
    sec.set_xlabel("f (Hz)")
    sec.xaxis.set_major_formatter(FuncFormatter(_hz_label))

    ax.set_ylim(0, max(1.0, mag.max()) * 1.1)
    ax.set_xlabel("Ω (rad/s)")
    ax.set_ylabel("|H(jΩ)|")
    ax.set_title(f"{title} — |H(jΩ)| (analog prototype)")
    ax.legend(loc="upper right")
    fig.tight_layout()
    return fig


# ------------------------------------------------- 3. z-plane: poles, zeros, ROC

def z_plane(b, a, title):
    z, p, _ = signal.tf2zpk(b, a)
    rmax = np.abs(p).max()
    lim = 1.6

    fig, ax = plt.subplots(figsize=(7, 7))

    # ROC (causal system): |z| > max|p| → everything except a disk of radius rmax
    ax.set_facecolor((0.17, 0.63, 0.17, 0.12))
    ax.add_patch(Circle((0, 0), rmax, color="white", zorder=0))
    ax.add_patch(Circle((0, 0), rmax, fill=False, ls=":", color="tab:green", lw=1.2,
                        label=f"ROC: |z| > {rmax:.3f}"))
    ax.add_patch(Circle((0, 0), 1, fill=False, color="k", lw=1.6, label="|z| = 1 (Fourier)"))

    ax.axhline(0, color="k", lw=0.6)
    ax.axvline(0, color="k", lw=0.6)

    _plot_roots(ax, p, 1.0, "x", "crimson", "poles")
    _plot_roots(ax, z, 1.0, "o", "royalblue", "zeros")

    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_xlabel("Re{z}")
    ax.set_ylabel("Im{z}")
    ax.set_title(f"{title} — z-plane")
    ax.legend(loc="upper right", fontsize=8)
    return fig


# ------------------------------------------- 4. Digital frequency response H(e^jω)

def _pi_label(x, _pos=None):
    n = round(x / (np.pi / 4))  # multiples of π/4
    labels = {0: "0", 4: "π", -4: "−π", 2: "π/2", -2: "−π/2"}
    if n in labels:
        return labels[n]
    num, den = (n, 4) if n % 2 else (n // 2, 2)
    sign = "−" if num < 0 else ""
    num = abs(num)
    return f"{sign}{'' if num == 1 else num}π/{den}"


def freq_digital(b, a, fs, f0, title):
    w = np.linspace(-np.pi, np.pi, 8193)
    _, h = signal.freqz(b, a, worN=w)
    mag = np.abs(h)
    w0 = 2 * np.pi * f0 / fs
    _, h0 = signal.freqz(b, a, worN=[w0])
    m0 = np.abs(h0[0])

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(w, mag)
    ax.axvline(w0, color="crimson", ls="--", lw=1, label=f"±ω0 = ±{w0 / np.pi:.4f}·π")
    ax.axvline(-w0, color="crimson", ls="--", lw=1)
    ax.scatter([-w0, w0], [m0, m0], color="crimson", zorder=5)
    ax.annotate(f"|H(e^jω0)| = {m0:.3f}", (w0, m0), xytext=(10, 6), textcoords="offset points",
                color="crimson")

    ax.set_xlim(-np.pi, np.pi)
    ax.xaxis.set_major_locator(FixedLocator(np.arange(-4, 5) * np.pi / 4))
    ax.xaxis.set_major_formatter(FuncFormatter(_pi_label))
    ax.grid(True, alpha=0.3)

    sec = ax.secondary_xaxis("top", functions=(lambda x: x * fs / (2 * np.pi),
                                               lambda f: f * 2 * np.pi / fs))
    sec.set_xlabel("f (Hz)")
    sec.xaxis.set_major_formatter(FuncFormatter(_hz_label))

    ax.set_ylim(0, max(1.0, mag.max()) * 1.1)
    ax.set_xlabel("ω (rad/sample)")
    ax.set_ylabel("|H(e^jω)|")
    ax.set_title(f"{title} — |H(e^jω)|")
    ax.legend(loc="upper right")
    fig.tight_layout()
    return fig


# ------------------------------------------------------- 5. DAW-style view

def daw_view(b, a, fs, f0, q, label):
    bg, panel, grid, fg, curve = "#1b1d21", "#24272c", "#363a41", "#9aa0a8", "#ff9f1a"

    f = np.logspace(np.log10(20), np.log10(min(20_000, fs / 2 * 0.999)), 1500)
    _, h = signal.freqz(b, a, worN=f, fs=fs)
    mag = _db(h)
    _, h0 = signal.freqz(b, a, worN=[f0], fs=fs)

    fig, ax = plt.subplots(figsize=(11, 4.5))
    fig.patch.set_facecolor(bg)
    ax.set_facecolor(panel)

    ax.fill_between(f, mag, 0, color=curve, alpha=0.15, lw=0)
    ax.semilogx(f, mag, color=curve, lw=2.2)
    ax.axhline(0, color=fg, lw=0.8, alpha=0.6)
    ax.scatter([f0], [_db(h0)[0]], s=70, color=curve, edgecolor="white", zorder=5)
    ax.annotate(f"{_hz_label(f0)}Hz · Q {q:.2f}", (f0, _db(h0)[0]), xytext=(10, 10),
                textcoords="offset points", color="white", fontsize=9)

    ax.set_xlim(20, 20_000)
    ax.set_ylim(-30, 15)
    ax.xaxis.set_major_locator(FixedLocator([20, 50, 100, 200, 500, 1000, 2000, 5000, 10_000, 20_000]))
    ax.xaxis.set_major_formatter(FuncFormatter(_hz_label))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_yticks(np.arange(-30, 16, 6))
    ax.grid(True, which="major", color=grid, lw=0.8)
    ax.grid(True, which="minor", axis="x", color=grid, lw=0.4, alpha=0.5)

    ax.tick_params(colors=fg, labelsize=8)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_ylabel("dB", color=fg)
    ax.text(0.01, 0.95, label, transform=ax.transAxes, color="white", fontsize=11,
            fontweight="bold", va="top")
    ax.text(0.99, 0.95, f"fs {fs / 1000:g} kHz", transform=ax.transAxes, color=fg,
            fontsize=8, va="top", ha="right")
    fig.tight_layout()
    return fig