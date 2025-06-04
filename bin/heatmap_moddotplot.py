import argparse
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Ensure text in PDF output is stored as text, not as shapes
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = False  # Disable LaTeX rendering

def main():
    parser = argparse.ArgumentParser(description="Heatmap of sequence identity values (symmetric)")
    parser.add_argument("input_file", help="Input file with 7 columns")
    parser.add_argument("output_file", help="Output filename")
    parser.add_argument("-pdf", action='store_true', help="Output in PDF format")
    parser.add_argument("--palette-colors", nargs='+', required=True, help="Hex colors for heatmap scale")
    parser.add_argument("--breakpoints", nargs='+', type=float, required=True, help="Breakpoints for color bins")
    parser.add_argument("--upper-triangle", action='store_true', help="Only show upper triangle")
    parser.add_argument("--xlim", nargs=2, type=int, help="X-axis limits")
    parser.add_argument("--ylim", nargs=2, type=int, help="Y-axis limits")
    parser.add_argument("--width", type=float, default=8.0, help="Width of output figure (inches)")
    parser.add_argument("--height", type=float, default=8.0, help="Height of output figure (inches)")
    parser.add_argument("-xlab", type=str, default="Position", help="Label for x-axis")
    parser.add_argument("-ylab", type=str, default="Position", help="Label for y-axis")
    parser.add_argument("--hide-scale-bar", action='store_true', help="Hide the identity color scale bar")
    parser.add_argument("--title", type=str, default="", help="Title of the whole plot")

    args = parser.parse_args()

    if len(args.palette_colors) != len(args.breakpoints) - 1:
        raise ValueError("Number of palette colors must be one less than number of breakpoints.")

    # Read data
    df = pd.read_csv(args.input_file, sep='\t', comment='#', header=None,
                     names=["query_name", "query_start", "query_end", "reference_name", "reference_start", "reference_end", "perID_by_events"])
    df['query_mid'] = (df['query_start'] + df['query_end']) // 2
    df['ref_mid'] = (df['reference_start'] + df['reference_end']) // 2

    if args.upper_triangle:
        df = df[df['query_mid'] <= df['ref_mid']]
    else:
        # Mirror the data to show the full square
        mirrored_df = df[df['query_mid'] != df['ref_mid']].copy()
        mirrored_df = mirrored_df.rename(columns={
            'query_mid': 'ref_mid',
            'ref_mid': 'query_mid'
        })
        df = pd.concat([df, mirrored_df], ignore_index=True)

    # Prepare colormap
    cmap = mcolors.ListedColormap(args.palette_colors)
    norm = mcolors.BoundaryNorm(args.breakpoints, cmap.N)

    # Plot
    fig, ax = plt.subplots(figsize=(args.width, args.height))
    sc = ax.scatter(df['query_mid'], df['ref_mid'], c=df['perID_by_events'], cmap=cmap, norm=norm, s=4, marker='s')

    # Axis config
    ax.set_xlabel(args.xlab)
    ax.set_ylabel(args.ylab)
    ax.set_aspect('equal', adjustable='box')

    if args.xlim:
        ax.set_xlim(args.xlim)
    else:
        ax.set_xlim(df['query_mid'].min(), df['query_mid'].max())

    if args.ylim:
        ax.set_ylim(args.ylim)
    else:
        ax.set_ylim(df['ref_mid'].min(), df['ref_mid'].max())

    # Move axes to origin
    ax.spines['left'].set_position('zero')
    ax.spines['bottom'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    # Add colorbar unless hidden
    if not args.hide_scale_bar:
        cbar = plt.colorbar(sc, ax=ax, boundaries=args.breakpoints, ticks=args.breakpoints)
        cbar.set_label("Identity (%)")

    # Add plot title
    if args.title:
        plt.title(args.title)

    plt.tight_layout()
    fmt = 'pdf' if args.pdf else args.output_file.split('.')[-1]
    plt.savefig(args.output_file, format=fmt)
    plt.close()

if __name__ == "__main__":
    main()

