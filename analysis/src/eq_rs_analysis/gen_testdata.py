from pathlib import Path

import numpy as np
from scipy.io import wavfile
from scipy.signal import chirp

OUT = Path(__file__).resolve().parents[3] / "testdata"
DUR = 3.0  # duration in s
AMP = 0.5 


def save(name: str, x: np.ndarray, fs: int) -> None:
    """Save signal in float32 and int16."""
    OUT.mkdir(exist_ok=True)
    x = x.astype(np.float32)
    wavfile.write(OUT / f"{name}_{fs}_f32.wav", fs, x)
    wavfile.write(OUT / f"{name}_{fs}_i16.wav", fs, np.round(x * 32767).astype(np.int16))


for fs in (44_100, 48_000):
    t = np.arange(int(DUR * fs)) / fs

    for f in (30, 1_000, 10_000):
        save(f"sine_{f}hz", AMP * np.sin(2 * np.pi * f * t), fs)

    rng = np.random.default_rng(0)
    save("noise", AMP * rng.uniform(-1, 1, t.size), fs)

    save("sweep", AMP * chirp(t, f0=20, t1=DUR, f1=20_000, method="logarithmic"), fs)

    impulse = np.zeros(fs)  # 1 s
    impulse[0] = 1.0
    save("impulse", impulse, fs)

    # Stereo: 1 kHz on L, 10 kHz on R
    stereo = AMP * np.column_stack([np.sin(2 * np.pi * 1_000 * t), np.sin(2 * np.pi * 10_000 * t)])
    save("stereo_1k_10k", stereo, fs)