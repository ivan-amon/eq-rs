# eq-rs · analysis

The math behind the filters in `eq-rs`.

Every filter in the Rust crate is a biquad from Robert Bristow-Johnson's
[Audio EQ Cookbook](https://www.w3.org/TR/audio-eq-cookbook/). This package derives
each filter from its analog prototype, checks that the cookbook coefficients are correct,
and plots the filter in the s-plane, the z-plane and the frequency domain.

Nothing here is used by the Rust code at runtime. It is only for reading and checking the math.

> [!IMPORTANT]
> **Notation: $\omega$ for continuous time, $\Omega$ for discrete time.** This repo uses lowercase
> $\omega$ for analog frequency and uppercase $\Omega$ for digital frequency.
> Many DSP references (Oppenheim & Schafer's *Discrete-Time Signal Processing*, Proakis & Manolakis,
> the scipy docs, the RBJ Cookbook) swap them: $\Omega$ for analog and $\omega$ for digital.
>
> | Symbol | Meaning | Units |
> |---|---|---|
> | $\omega$ | analog (continuous-time) angular frequency, as in $H(j\omega)$ | rad/s |
> | $\omega_c$ | cutoff of the analog prototype (after prewarping) | rad/s |
> | $\Omega$ | digital (discrete-time) angular frequency, as in $H(e^{j\Omega})$ | rad/sample |
> | $\Omega_0 = \omega_0 T_s$ | digital cutoff / center frequency, with $\omega_0 = 2\pi f_0$ and $T_s = 1/f_s$ | rad/sample |
> | $f_s$ | sample rate | Hz |
>
> The RBJ Cookbook calls the digital frequency `w0`. The code keeps the name `w0` so it can be
> compared line by line with the cookbook, but in this notation it means $\Omega_0$.

## Quick start

You need [uv](https://docs.astral.sh/uv/). It installs Python 3.14 and the dependencies for you.

```sh
cd analysis
uv run lpf              # every plot for the default LPF (f0 = 1 kHz, Q = 0.707, fs = 48 kHz)
uv run lpf --f0 300 --q 4
uv run lpf -d           # digital plots only
uv run lpf --save figs/ # write PNGs instead of opening windows
```

## Options

| Option | Default | Meaning |
|---|---|---|
| `--f0` | `1000` | cutoff / center frequency (Hz) |
| `--q` | `1/√2` | quality factor |
| `--fs` | `48000` | sample rate (Hz) |
| `--gain` | `6` | gain in dB (peaking and shelving filters only) |
| `--save DIR` | | save PNGs to `DIR` instead of showing windows |

To show only some of the plots, use the view flags below. They can be combined
(for example `-a --daw`). With no view flag, every plot is shown.

| Flag | Plots |
|---|---|
| `-a`, `--analog` | s-plane, analog frequency response |
| `-d`, `--digital` | z-plane, digital frequency response |
| `--daw` | DAW-style view |

## What each plot shows

1. **s-plane**: the poles and zeros of the analog prototype $H(s)$, with the region of convergence.
2. **Analog frequency response**: $H(j\omega)$, which is the Laplace transform evaluated on the imaginary axis.
   It covers $\omega \in [-\pi f_s, \pi f_s]$, the same span as the digital plot.
3. **z-plane**: the poles and zeros of $H(z)$, with the unit circle and the region of convergence.
4. **Digital frequency response**: $H(e^{j\Omega})$ for $\Omega \in [-\pi, \pi]$, which is $H(z)$ evaluated on the unit circle.
5. **DAW view**: magnitude in dB on a log axis from 20 Hz to 20 kHz, styled like the EQ in an audio plugin.

Each run also prints the coefficients `b` and `a`, the prewarped $\omega_c$ and the z-plane poles.

## From analog prototype to biquad

The same steps apply to every filter:

1. **Analog prototype.** Start from a second-order $H(s)$. For example, the LPF is

   $$H(s) = \frac{\omega_c^2}{s^2 + \frac{\omega_c}{Q}s + \omega_c^2}$$

   With $Q = 1/\sqrt{2}$, this is a 2nd-order Butterworth filter.

2. **Prewarping.** The bilinear transform compresses the frequency axis, so the analog cutoff
   is chosen to land exactly on $f_0$ after the transform:

   $$\omega_c = \frac{2}{T_s} \tan\left(\frac{\Omega_0}{2}\right) = \frac{2}{T_s} \tan\left(\frac{\omega_0 T_s}{2}\right)$$

3. **Bilinear transform.** Substitute

   $$s = \frac{2}{T_s} \cdot \frac{1 - z^{-1}}{1 + z^{-1}}$$

   This maps the left half of the s-plane into the unit circle, so a stable analog filter
   gives a stable digital filter.

4. **Cookbook coefficients.** The result is a biquad, normalized so that $a_0 = 1$:

   $$H(z) = \frac{b_0 + b_1 z^{-1} + b_2 z^{-2}}{1 + a_1 z^{-1} + a_2 z^{-2}}$$

On every run, an assertion checks that `scipy.signal.bilinear(prototype)` equals the cookbook
coefficients. If someone mistypes a formula, the run fails before anything is plotted.

## Conventions

These match `scipy.signal` and are the conventions the Rust code should follow:

- Analog `num`, `den`: polynomials in $s$, highest power first (as in `scipy.signal.freqs`).
- Digital `b`, `a`: coefficients in $z^{-1}$, with `a[0] = 1` (as in `scipy.signal.freqz`).
- Filters are causal, so the region of convergence is to the right of the rightmost pole in $s$
  and outside the outermost pole in $z$.
- Frequency symbols: $\omega$ is analog (rad/s), $\Omega$ is digital (rad/sample), as described at the top of this README.

## Layout

```
src/eq_rs_analysis/
├── rbj.py    # cookbook helpers: Ω0, α, A, prewarping, a0 normalization
├── cli.py    # shared runner: CLI arguments, bilinear check, choosing plots
├── plots.py  # generic plots for any H(s) / H(z)
└── lpf.py    # low-pass filter (one module per filter)
```

## Adding a filter

1. Create `src/eq_rs_analysis/<name>.py` with two functions, using `lpf.py` as the template:
   - `digital(f0, q, fs, gain_db) -> (b, a)`: the cookbook coefficients, normalized with `rbj.normalize`.
   - `analog(f0, q, fs, gain_db) -> (num, den, ωc)`: the prototype, with $\omega_c$ from `rbj.prewarp`.
2. Define a `Filter(...)` with them and call `run(FILTER)` from `main()`.
   Pass `has_gain=True` for peaking and shelving filters, which adds `--gain`.
3. Register the command in `pyproject.toml`:

   ```toml
   [project.scripts]
   hpf = "eq_rs_analysis.hpf:main"
   ```

## References

- [`docs/biquads.md`](../docs/biquads.md): the formulas of every filter, in this repo's notation.
- R. Bristow-Johnson, [*Audio EQ Cookbook*](https://www.w3.org/TR/audio-eq-cookbook/).
- [`scipy.signal`](https://docs.scipy.org/doc/scipy/reference/signal.html): `bilinear`, `freqs`, `freqz`, `tf2zpk`.
