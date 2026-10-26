#!/usr/bin/env python3
"""
plot_displacement_pdf.py

Takes several *_traj.dat trajectory files and plots the probability 
density function (PDF) of node displacements relative to their time-averaged 
positions: P(x - <x>) and P(y - <y>).

Usage:
    python3 plot_displacement_pdf.py traj_files/*_traj.dat --out pdf_x.png
    
    # Plot the y-displacements with a log scale on the y-axis:
    python3 plot_displacement_pdf.py traj_files/*_traj.dat --axis y --logy --out pdf_y.png
"""

import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from msd_utils import parse_trajectory, unwrap_positions
# Import the new function we just added
from msd_utils_2 import compute_displacement_deviations

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("traj_files", nargs="+", help="*_traj.dat files, one per z")
    ap.add_argument("--out", default="displacement_pdf.png", help="output plot path")
    ap.add_argument("--axis", choices=["x", "y", "both"], default="x",
                    help="Which component to histogram (default: x)")
    ap.add_argument("--bins", type=int, default=100, help="Number of histogram bins")
    ap.add_argument("--logy", action="store_true", 
                    help="Log-scale the PDF axis (useful for spotting exponential tails)")
    args = ap.parse_args(argv)

    results = []
    
    for path in args.traj_files:
        N, contacts, L, strain, bonds, steps, frames = parse_trajectory(path)
        z = 2.0 * contacts / N
        
        unwrapped = unwrap_positions(frames, L)
        
        # Calculate x - <x> and y - <y>
        dx, dy = compute_displacement_deviations(unwrapped)
        
        print(f"Processed {path}: z={z:.4g}, N={N}, frames={len(frames)}")
        results.append((z, dx, dy))

    # Sort by z-coordination to keep the legend organized
    results.sort(key=lambda r: r[0])

    fig, ax = plt.subplots(figsize=(7, 5))
    cmap = plt.get_cmap("viridis")

    for k, (z, dx, dy) in enumerate(results):
        color = cmap(k / max(1, len(results) - 1))
        
        data_to_plot = []
        if args.axis in ["x", "both"]:
            data_to_plot.append((dx, f"z={z:.3f} (x)"))
        if args.axis in ["y", "both"]:
            data_to_plot.append((dy, f"z={z:.3f} (y)"))
            
        for data, label in data_to_plot:
            # density=True ensures the area under the histogram integrates to 1 (PDF)
            # histtype='step' draws outlines rather than filled bars so overlapping curves are visible
            ax.hist(data, bins=args.bins, density=True, histtype='step', 
                    color=color, lw=1.5, label=label, alpha=0.8)

    ax.set_xlabel(r"Displacement ($x - \langle x \rangle$ or $y - \langle y \rangle$)")
    ax.set_ylabel("Probability Density")
    ax.set_title(f"Displacement PDF across different z")
    
    if args.logy:
        ax.set_yscale("log")
        # Prevent log(0) errors on the plot by setting a sensible lower limit
        ax.set_ylim(bottom=1e-5) 
        
    ax.legend(fontsize=8, loc='upper right', bbox_to_anchor=(1.25, 1.0))
    fig.tight_layout()
    fig.savefig(args.out, dpi=150, bbox_inches='tight')
    print(f"Wrote {args.out}")

if __name__ == "__main__":
    main()