#!/usr/bin/env python3
"""
make_network_movie.py

Turns the trajectory file produced by the patched rigid_bar_network_2d.c
(runNVE now writes "<network>_traj.dat" alongside the usual
"<network>_nve.dat" diagnostics file) into an mp4 or gif movie of the
network's motion.

Trajectory file format (written by dumpFrame() in the C code):

    # TRAJ N=<N> contacts=<Nc> L=<L> strain=<strain>
    # BONDS
    i_0   j_0
    i_1   j_1
    ...                      (Nc lines, bond topology -- fixed for the run)
    # FRAMES
    FRAME <step_0>
    rx_0  ry_0
    rx_1  ry_1
    ...                      (N lines, positions in reduced [0,1] units)
    FRAME <step_1>
    ...

Usage:
    python make_network_movie.py floppy_network_N1600_z2_800_s00000_traj.dat \
        --out network_movie.mp4 --fps 15

If ffmpeg isn't available, pass --gif to fall back to a Pillow-written GIF.
"""

import argparse
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter


def parse_trajectory(path):
    """Parse the # TRAJ / # BONDS / # FRAMES file into (N, contacts, L,
    strain, bonds[Nc,2] int, steps[list], frames[list of (N,2) arrays])."""
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
                # "# TRAJ N=1600 contacts=2240 L=40 strain=0"
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
                    frames.append(np.array(cur_frame, dtype=float))
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
        frames.append(np.array(cur_frame, dtype=float))

    bonds = np.array(bonds, dtype=int)
    if N is None or contacts is None:
        raise ValueError("Could not find '# TRAJ ...' header line -- "
                          "is this really a *_traj.dat file from the "
                          "patched rigid_bar_network_2d.c?")
    return N, contacts, L, strain, bonds, steps, frames


def minimum_image_segments(pos, bonds, strain):
    """Build the (x0,y0)-(x1,y1) line segments for every bond, applying
    the SAME Lees-Edwards minimum-image convention used in the C code's
    readNetwork(), so bonds that wrap around the periodic box are drawn
    as the short way around instead of a spurious line across the cell.

    pos: (N,2) array in reduced [0,1] units.
    Returns an (Nc, 2, 2) array suitable for LineCollection.
    """
    le = strain % 1.0
    if le < -0.5:
        le += 1.0
    elif le > 0.5:
        le -= 1.0

    i = bonds[:, 0]
    j = bonds[:, 1]
    dx = pos[j, 0] - pos[i, 0]
    dy = pos[j, 1] - pos[i, 1]

    wrap_hi = dy >= 0.5
    wrap_lo = dy < -0.5
    dy[wrap_hi] -= 1.0
    dx[wrap_hi] -= le
    dy[wrap_lo] += 1.0
    dx[wrap_lo] += le

    dx[dx >= 0.5] -= 1.0
    dx[dx < -0.5] += 1.0

    x0, y0 = pos[i, 0], pos[i, 1]
    x1, y1 = x0 + dx, y0 + dy
    segs = np.stack([np.stack([x0, y0], axis=1),
                      np.stack([x1, y1], axis=1)], axis=1)
    return segs


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("traj_file", help="*_traj.dat file from runNVE")
    ap.add_argument("--out", default="network_movie.mp4",
                     help="output movie path (default: network_movie.mp4)")
    ap.add_argument("--fps", type=int, default=15)
    ap.add_argument("--gif", action="store_true",
                     help="force GIF output via Pillow instead of ffmpeg")
    ap.add_argument("--dpi", type=int, default=140)
    ap.add_argument("--node-size", type=float, default=4.0)
    args = ap.parse_args(argv)

    N, contacts, L, strain, bonds, steps, frames = parse_trajectory(args.traj_file)
    n_frames = len(frames)
    print(f"Parsed {args.traj_file}: N={N}  contacts={contacts}  "
          f"L={L:.4g}  strain={strain:.4g}  frames={n_frames}")
    if n_frames < 2:
        sys.exit("Need at least 2 frames to make a movie -- rerun the C "
                  "code with more steps or a smaller out_every.")

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_facecolor("#0b1020")
    fig.patch.set_facecolor("#0b1020")

    lc = LineCollection([], linewidths=0.6, colors="#4fd1c5", alpha=0.85)
    ax.add_collection(lc)
    scat = ax.scatter([], [], s=args.node_size, c="#f6f7fb",
                       linewidths=0, zorder=3)
    title = ax.text(0.5, 1.02, "", transform=ax.transAxes, ha="center",
                     va="bottom", color="#e5e7eb", fontsize=10)

    def init():
        lc.set_segments([])
        scat.set_offsets(np.zeros((N, 2)))
        title.set_text("")
        return lc, scat, title

    def update(frame_idx):
        pos = frames[frame_idx]
        segs = minimum_image_segments(pos, bonds, strain)
        lc.set_segments(segs)
        scat.set_offsets(pos)
        title.set_text(f"step {steps[frame_idx]}")
        return lc, scat, title

    anim = FuncAnimation(fig, update, frames=n_frames, init_func=init,
                          blit=True, interval=1000 / args.fps)

    out_path = args.out
    use_gif = args.gif or out_path.lower().endswith(".gif")
    if use_gif:
        if not out_path.lower().endswith(".gif"):
            out_path += ".gif"
        writer = PillowWriter(fps=args.fps)
    else:
        writer = FFMpegWriter(fps=args.fps, bitrate=2400)

    print(f"Writing {out_path} ({n_frames} frames @ {args.fps} fps)...")
    anim.save(out_path, writer=writer, dpi=args.dpi)
    plt.close(fig)
    print("Done:", out_path)


if __name__ == "__main__":
    main()
