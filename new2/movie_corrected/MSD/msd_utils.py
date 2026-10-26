"""
msd_utils.py

Shared helpers for computing mean-squared displacement (MSD) from the
*_traj.dat trajectory files written by the patched rigid_bar_network_2d.c
(same format for runNVE and runDriftTest movie dumps):

    # TRAJ N=<N> contacts=<Nc> L=<L> strain=<strain>
    # BONDS
    i_0   j_0
    ...
    # FRAMES
    FRAME <step_0>
    rx_0  ry_0
    ...
    FRAME <step_1>
    ...

Positions are reduced coordinates in [0,1]; multiply by L for actual units.
"""

import numpy as np


def parse_trajectory(path):
    """Returns (N, contacts, L, strain, bonds[Nc,2] int, steps[list int],
    frames: (n_frames, N, 2) float array in reduced [0,1] units)."""
    N = contacts = None
    L = strain = None
    bonds = []
    steps = []
    frames = []
    cur_frame = None
    section = None

    with open(path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("# TRAJ"):
                parts = dict(tok.split("=") for tok in line.split()[2:])
                N = int(parts["N"])
                contacts = int(parts["contacts"])
                L = float(parts["L"])
                strain = float(parts["strain"])
                continue
            if line.startswith("# BONDS"):
                section = "bonds"
                continue
            if line.startswith("# FRAMES"):
                section = "frames"
                continue
            if line.startswith("FRAME"):
                if cur_frame is not None:
                    frames.append(cur_frame)
                steps.append(int(line.split()[1]))
                cur_frame = []
                continue
            if section == "bonds":
                i, j = line.split()
                bonds.append((int(i), int(j)))
            elif section == "frames":
                rx, ry = line.split()
                cur_frame.append((float(rx), float(ry)))

    if cur_frame is not None:
        frames.append(cur_frame)

    if N is None or contacts is None:
        raise ValueError(f"{path}: no '# TRAJ ...' header found -- is this "
                          f"really a *_traj.dat file?")

    bonds = np.array(bonds, dtype=int)
    steps = np.array(steps, dtype=int)
    frames = np.array(frames, dtype=float)  # (n_frames, N, 2)
    return N, contacts, L, strain, bonds, steps, frames


def unwrap_positions(frames, L):
    """Frame-to-frame minimum-image unwrapping (NOT bond-based -- this is
    per-node, since MSD tracks individual nodes, not bonds).

    Assumes consecutive saved frames are close enough in time that no
    node moves more than half the box in one recorded interval; if your
    out_every is very large relative to how fast the network moves, this
    assumption breaks and displacements alias. Rerun with a smaller
    out_every / frame_every if MSD looks suspiciously small or noisy.

    frames: (n_frames, N, 2) reduced-unit positions.
    Returns unwrapped positions in ACTUAL units (already multiplied by L),
    shape (n_frames, N, 2).
    """
    n_frames, N, _ = frames.shape
    raw_disp = np.diff(frames, axis=0)          # (n_frames-1, N, 2), reduced units
    raw_disp[raw_disp > 0.5]  -= 1.0
    raw_disp[raw_disp < -0.5] += 1.0

    unwrapped = np.empty_like(frames)
    unwrapped[0] = frames[0]
    unwrapped[1:] = frames[0] + np.cumsum(raw_disp, axis=0)
    return unwrapped * L


def compute_msd(unwrapped, max_lag=None, subtract_com=True):
    """Time-averaged MSD over all valid time origins for each lag.

    unwrapped: (n_frames, N, 2) actual-unit unwrapped positions.
    subtract_com: if True, subtracts the network's center-of-mass motion
        at each frame before computing MSD, so overall rigid translation
        (which should be ~0 anyway if momentum is conserved, but isn't
        exactly zero at floating-point precision) doesn't contaminate the
        internal-rearrangement signal.
    max_lag: largest lag (in frames) to compute; defaults to n_frames // 2,
        the usual rule of thumb so each lag still averages over enough
        time origins to be statistically meaningful.

    Returns (lags: (max_lag,) int array of frame-lags, msd: (max_lag,) array).
    """
    n_frames = unwrapped.shape[0]
    if max_lag is None:
        max_lag = max(1, n_frames // 2)
    max_lag = min(max_lag, n_frames - 1)

    pos = unwrapped
    if subtract_com:
        pos = pos - pos.mean(axis=1, keepdims=True)

    lags = np.arange(1, max_lag + 1)
    msd = np.empty(max_lag, dtype=float)
    for k, lag in enumerate(lags):
        diff = pos[lag:] - pos[:-lag]              # (n_frames-lag, N, 2)
        sq = np.sum(diff**2, axis=-1)               # (n_frames-lag, N)
        msd[k] = sq.mean()
    return lags, msd
