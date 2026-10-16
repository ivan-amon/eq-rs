# Biquad filters

The formulas behind every filter in `eq-rs`. The Rust code implements them, and the Python
package in [`analysis/`](../analysis/) checks them and plots them.

## Notation

The whole repo uses lowercase $\omega$ for continuous time and uppercase $\Omega$ for discrete time.

| Symbol | Meaning | Units |
|---|---|---|
| $f_s$ | sample rate | Hz |
| $T_s = 1 / f_s$ | sampling period | s |
| $f_0$ | cutoff / center frequency | Hz |
| $\omega$ | analog angular frequency, as in $H(j\omega)$ | rad/s |
| $\omega_0 = 2\pi f_0$ | cutoff / center frequency | rad/s |
| $\omega_c$ | cutoff of the analog prototype, after prewarping | rad/s |
| $\Omega$ | digital angular frequency, as in $H(e^{j\Omega})$ | rad/sample |
| $\Omega_0 = \omega_0 T_s$ | digital cutoff / center frequency | rad/sample |
| $Q$ | quality factor | |
| $s = \sigma + j\omega$ | Laplace variable | |
| $z = r e^{j\Omega}$ | Z-transform variable | |

The two frequencies are related by $\Omega = \omega T_s$ for a sampled signal.

Many DSP references, including the Audio EQ Cookbook, swap the two symbols:
they write $\Omega$ for analog frequency and $\omega$ for digital frequency.
Where the cookbook writes `w0`, this repo writes $\Omega_0$. The Python code keeps the variable
name `w0` so it can be compared line by line with the cookbook.

## The biquad

Every filter is a second-order IIR filter. Its coefficients are normalized so that $a_0 = 1$:

$$H(z) = \frac{b_0 + b_1 z^{-1} + b_2 z^{-2}}{1 + a_1 z^{-1} + a_2 z^{-2}}$$

The formulas below give un-normalized coefficients. Divide all six by $a_0$ before use.

In the time domain this is the difference equation

$$y[n] = b_0 x[n] + b_1 x[n-1] + b_2 x[n-2] - a_1 y[n-1] - a_2 y[n-2]$$

## Design method

Each filter starts as a second-order analog prototype $H(s)$ and is converted to $H(z)$
with the bilinear transform.

1. **Prewarping.** The bilinear transform compresses the frequency axis. The analog cutoff is
   chosen so that it lands exactly on $\Omega_0$ after the transform:

   $$\omega_c = \frac{2}{T_s} \tan\left(\frac{\Omega_0}{2}\right)$$

   Since $\tan x > x$, $\omega_c > \omega_0$. The two are almost equal when $f_0 \ll f_s / 2$
   (1 kHz at 48 kHz gives 1001.4 Hz) and grow apart as $f_0$ approaches $f_s / 2$.

2. **Bilinear transform.** Substitute

   $$s = \frac{2}{T_s} \cdot \frac{1 - z^{-1}}{1 + z^{-1}}$$

   The left half of the s-plane maps to the inside of the unit circle, so a stable prototype
   always gives a stable biquad.

The cookbook writes the prototypes with a normalized cutoff ($\omega_c = 1$). With prewarping,
both forms give the same result, because

$$\frac{s}{\omega_c} = \frac{1}{\tan(\Omega_0 / 2)} \cdot \frac{1 - z^{-1}}{1 + z^{-1}}$$

## Shared parameters

$$\Omega_0 = \omega_0 T_s \qquad \alpha = \frac{\sin \Omega_0}{2Q}$$

**Stability.** The biquad is stable if and only if $\alpha > 0$. That means:

- $Q > 0$
- $0 < \Omega_0 < \pi$, that is, $0 < f_0 < f_s / 2$

## LPF (low-pass)

**Analog prototype:**

$$H(s) = \frac{\omega_c^2}{s^2 + \frac{\omega_c}{Q}s + \omega_c^2}$$

**Coefficients:**

$$
\begin{aligned}
b_0 &= \frac{1 - \cos \Omega_0}{2} & a_0 &= 1 + \alpha \\
b_1 &= 1 - \cos \Omega_0 & a_1 &= -2 \cos \Omega_0 \\
b_2 &= \frac{1 - \cos \Omega_0}{2} & a_2 &= 1 - \alpha
\end{aligned}
$$

**Properties:**

- Unity gain at DC: $H(e^{j0}) = 1$.
- Double zero at $z = -1$, so the gain at Nyquist is $0$.
- The gain at the cutoff equals $Q$: $|H(e^{j\Omega_0})| = Q$.
  With $Q = 1/\sqrt{2}$, the filter is a 2nd-order Butterworth and the cutoff is the $-3$ dB point.
- Slope of 12 dB/octave above the cutoff.

**Checked in:** [`analysis/src/eq_rs_analysis/lpf.py`](../analysis/src/eq_rs_analysis/lpf.py)
(`uv run lpf`).

## Source

The coefficients come from R. Bristow-Johnson,
[*Audio EQ Cookbook*](https://www.w3.org/TR/audio-eq-cookbook/), W3C Working Group Note,
8 June 2021. They are rewritten here in this repo's notation, with the analog prototypes
written with an explicit $\omega_c$.
