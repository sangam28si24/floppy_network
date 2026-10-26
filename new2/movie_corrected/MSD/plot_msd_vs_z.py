#!/usr/bin/env python3
"""
plot_msd_vs_z.py

Takes several *_traj.dat trajectory files -- one per network at a
different coordination number z -- reduces each one's MSD(lag) curve to
a single summary number, and plots that summary MSD against z.

z is read directly from each file's own "# TRAJ ... contacts=... N=..."
header (z = 2*contacts/N), not parsed from the filename, so it doesn't
matter how you name your files.

Two summary modes:
  --summary plateau (default): average MSD over the last `--plateau-frac`
      fraction of the lag range (e.g. 0.3 = last 30% of lags). Good for a
      long-time / saturation value.
  --summary lag: MSD evaluated at one specific lag (in steps), given by
      --lag-steps. Use this to compare all z values at the *same* physical
      lag time rather than each curve's own tail.

Usage:
    python3 plot_msd_vs_z.py traj_files/*_traj.dat --out msd_vs_z.png

    # compare all networks at exactly 2000 steps of lag:
    python3 plot_msd_vs_z.py traj_files/*_traj.dat \
        --summary lag --lag-steps 2000 --out msd_vs_z.png

    # also dump the individual MSD(lag) curves alongside the summary:
    python3 plot_msd_vs_z.py traj_files/*_traj.dat --out msd_vs_z.png \
        --curves-out msd_curves.png
"""

import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from msd_utils import parse_trajectory, unwrap_positions, compute_msd


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("traj_files", nargs="+", help="*_traj.dat files, one per z")
    ap.add_argument("--out", default="msd_vs_z.png", help="output plot path")
    ap.add_argument("--csv", default=None,
                     help="also write the summary table as CSV: z,msd_summary,N,contacts")
    ap.add_argument("--summary", choices=["plateau", "lag"], default="plateau")
    ap.add_argument("--plateau-frac", type=float, default=0.3,
                     help="fraction of the tail of each MSD(lag) curve to "
                          "average for --summary plateau (default 0.3)")
    ap.add_argument("--lag-steps", type=int, default=None,
                     help="required if --summary lag: the RATTLE-step lag "
                          "to evaluate MSD at for every file")
    ap.add_argument("--max-lag", type=int, default=None,
                     help="largest lag (in frames) to compute per file; "
                          "default n_frames//2")
    ap.add_argument("--curves-out", default=None,
                     help="optional second plot: all individual MSD(lag) "
                          "curves overlaid, colored by z")
    ap.add_argument("--logy", action="store_true",
                     help="log-scale the MSD axis on the summary plot")
    args = ap.parse_args(argv)

    if args.summary == "lag" and args.lag_steps is None:
        ap.error("--summary lag requires --lag-steps")

    results = []  # (z, summary_msd, N, contacts, lags_steps, msd_curve)
    for path in args.traj_files:
        N, contacts, L, strain, bonds, steps, frames = parse_trajectory(path)
        z = 2.0 * contacts / N
        unwrapped = unwrap_positions(frames, L)
        lags, msd = compute_msd(unwrapped, max_lag=args.max_lag)
        lag_steps = lags * (steps[1] - steps[0])

        if args.summary == "plateau":
            n_tail = max(1, int(len(msd) * args.plateau_frac))
            summary = msd.mean()
        else:
            idx = int(np.argmin(np.abs(lag_steps - args.lag_steps)))
            summary = msd[idx]
            if abs(lag_steps[idx] - args.lag_steps) > 0.5 * (steps[1] - steps[0]):
                print(f"  [warn] {path}: requested lag {args.lag_steps} steps, "
                      f"nearest available is {lag_steps[idx]} steps")

        print(f"{path}: N={N} contacts={contacts} z={z:.4g} "
              f"n_frames={len(frames)} -> summary MSD = {summary:.6g}")
        results.append((z, summary, N, contacts, lag_steps, msd))

    results.sort(key=lambda r: r[0])
    zs = np.array([r[0] for r in results])
    summaries = np.array([r[1] for r in results])

    if args.csv:
        with open(args.csv, "w") as f:
            f.write("z\tmsd_summary\tN\tcontacts\n")
            for z, s, N, c, _, _ in results:
                f.write(f"{z}\t{s}\t{N}\t{c}\n")
        print("Wrote", args.csv)

    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(zs, summaries, "o-", lw=1.5, ms=6)
    ax.axvline(4.0, color="gray", ls="--", lw=1, label=r"$z_c=4$ (isostatic)")
    if args.logy:
        ax.set_yscale("log")
    ax.set_xlabel("z")
    label = (f"MSD"
             if args.summary == "plateau" else f"MSD @ lag={args.lag_steps} steps")
    ax.set_ylabel(label)
    ax.set_title("MSD vs z")
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print("Wrote", args.out)

    if args.curves_out:
        fig2, ax2 = plt.subplots(figsize=(6.5, 4.5))
        cmap = plt.get_cmap("viridis")
        for k, (z, s, N, c, lag_steps, msd) in enumerate(results):
            color = cmap(k / max(1, len(results) - 1))
            ax2.plot(lag_steps, msd, color=color, lw=1.3, label=f"z={z:.3f}")
        ax2.set_xlabel("lag (steps)")
        ax2.set_ylabel("MSD")
        ax2.set_title("MSD(lag) per z")
        ax2.legend(fontsize=8)
        fig2.tight_layout()
        fig2.savefig(args.curves_out, dpi=150)
        print("Wrote", args.curves_out)


if __name__ == "__main__":
    main()
