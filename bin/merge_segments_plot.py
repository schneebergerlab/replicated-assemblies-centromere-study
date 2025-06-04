import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import os

import matplotlib

# Ensure text in PDF output is stored as text, not as shapes
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = False  # Disable LaTeX rendering

# Argument parser
parser = argparse.ArgumentParser(description="Merge aligned segments and plot them.")
parser.add_argument("-i", "--input", required=True, help="Input file with alignment segments.")
parser.add_argument("-o", "--output", required=True, help="Output file for merged segment data.")
parser.add_argument("-pdf", action="store_true", help="Save plot as PDF instead of PNG.")
parser.add_argument("--highlight_reverse", action="store_true",
                    help="Highlight segments with negative slope (reverse alignment) in red.")
parser.add_argument("--xlab", default="Ref", help="Label for the x-axis.")
parser.add_argument("--ylab", default="Qry", help="Label for the y-axis.")
parser.add_argument("--flip", action="store_true",
                    help="Flip x and y axes in the plot.")
parser.add_argument("--xstart", type=int, default=0,
                    help="Starting coordinate for x-axis (default: 0).")
parser.add_argument("--ystart", type=int, default=0,
                    help="Starting coordinate for y-axis (default: 0).")
parser.add_argument("--label_internal_ends", action="store_true",
                    help="Label internal ends of the two blue-highlighted segments.")
args = parser.parse_args()

# Read input
data = pd.read_csv(args.input, sep="\t", header=0)

# Compute segment positions
data["x1"] = data["Pos_Seq1"]
data["x2"] = data["x1"] + data["Matched_Word"].apply(len) - 1
data["y1"] = np.where(data["Strand"] == "+",
                      data["Pos_Seq2"],
                      data["Pos_Seq2"] + data["Matched_Word"].apply(len) - 1)
data["y2"] = np.where(data["Strand"] == "+",
                      data["y1"] + data["Matched_Word"].apply(len) - 1,
                      data["y1"] - data["Matched_Word"].apply(len) + 1)

# Check if a point is on a segment
def point_on_segment(x, y, x1, y1, x2, y2):
    dx1, dy1 = x2 - x1, y2 - y1
    dx2, dy2 = x - x1, y - y1

    if dx1 * dy2 != dx2 * dy1:
        return False

    dot = dx1 * dx2 + dy1 * dy2
    if dot < 0:
        return False

    squared_len = dx1 * dx1 + dy1 * dy1
    return dot <= squared_len

# Merge segments
def merge_segments(df):
    merged = []
    for i in range(len(df)):
        x1, y1, x2, y2 = df.loc[i, ["x1", "y1", "x2", "y2"]]
        slope_curr = (y2 - y1) / (x2 - x1) if x2 != x1 else float("inf")
        merged_flag = False
        for seg in merged:
            mx1, my1, mx2, my2 = seg["x1"], seg["y1"], seg["x2"], seg["y2"]
            slope_seg = (my2 - my1) / (mx2 - mx1) if mx2 != mx1 else float("inf")
            if math.isclose(slope_curr, slope_seg, abs_tol=1e-6):
                if point_on_segment(x1, y1, mx1, my1, mx2, my2) or point_on_segment(x2, y2, mx1, my1, mx2, my2):
                    seg["x1"] = min(mx1, x1)
                    seg["y1"] = min(my1, y1) if slope_seg >= 0 else max(my1, y1)
                    seg["x2"] = max(mx2, x2)
                    seg["y2"] = max(my2, y2) if slope_seg >= 0 else min(my2, y2)
                    merged_flag = True
                    break
        if not merged_flag:
            merged.append({"x1": x1, "y1": y1, "x2": x2, "y2": y2})
    return pd.DataFrame(merged)

# Perform merging
merged_df = merge_segments(data)
merged_df["x1"] = merged_df["x1"] + args.xstart
merged_df["x2"] = merged_df["x2"] + args.xstart
merged_df["y1"] = merged_df["y1"] + args.ystart
merged_df["y2"] = merged_df["y2"] + args.ystart
merged_df["length"] = np.round(np.sqrt((merged_df["x2"] - merged_df["x1"])**2 + (merged_df["y2"] - merged_df["y1"])**2), 1)

# Save merged output
# Rename columns for clarity in output
output_df = merged_df.rename(columns={
    "x1": "Ref_start",
    "x2": "Ref_end",
    "y1": "Qry_start",
    "y2": "Qry_end",
    "length": "Length"
})
output_df.to_csv(args.output, sep="\t", index=False)
#merged_df.to_csv(args.output, sep="\t", index=False)

# Plot
# Find the bottom-left segment (based on x1 + y1)
lower_left_idx = merged_df.assign(key=merged_df["x1"] + merged_df["y1"]).sort_values("key").index[0]

# Find the top-right segment (based on x2 + y2)
upper_right_idx = merged_df.assign(key=merged_df["x2"] + merged_df["y2"]).sort_values("key", ascending=False).index[0]

# Set colors
# Assign colors based on slope and special positions
colors = []
for i, row in merged_df.iterrows():
    if i == lower_left_idx or i == upper_right_idx:
        colors.append("blue")  # Special endpoints
    else:
        slope = (row["y2"] - row["y1"]) / (row["x2"] - row["x1"]) if row["x2"] != row["x1"] else float("inf")
        if slope < 0 and args.highlight_reverse:
            colors.append("red")  # Highlight reverse alignment
        else:
            colors.append("darkgrey")  # Default color

#colors = ["blue" if i in [lower_left_idx, upper_right_idx] else "darkgrey" for i in merged_df.index]

#longest_idx = merged_df["length"].nlargest(2).index
#colors = ["blue" if i in longest_idx else "darkgrey" for i in merged_df.index]

plt.figure(figsize=(10, 8))
#for i, row in merged_df.iterrows():
#    plt.plot([row["x1"], row["x2"]], [row["y1"], row["y2"]], color=colors[i])
for i, row in merged_df.iterrows():
    if args.flip:
        plt.plot([row["y1"], row["y2"]], [row["x1"], row["x2"]], color=colors[i])
    else:
        plt.plot([row["x1"], row["x2"]], [row["y1"], row["y2"]], color=colors[i])
if args.flip:
    plt.xlabel(args.ylab)
    plt.ylabel(args.xlab)
else:
    plt.xlabel(args.xlab)
    plt.ylabel(args.ylab)

#plt.title("Merged Alignment Segments")
#plt.grid(True)

if args.label_internal_ends:
    # Get the two special blue-highlighted segments
    seg1 = merged_df.loc[lower_left_idx]  # Bottom-left segment
    seg2 = merged_df.loc[upper_right_idx]  # Top-right segment

    # For bottom-left segment, internal end = right-up corner
    if seg1["x1"] + seg1["y1"] < seg1["x2"] + seg1["y2"]:
        x_int1, y_int1 = seg1["x2"], seg1["y2"]
    else:
        x_int1, y_int1 = seg1["x1"], seg1["y1"]

    # For top-right segment, internal end = left-lower corner
    if seg2["x1"] + seg2["y1"] < seg2["x2"] + seg2["y2"]:
        x_int2, y_int2 = seg2["x1"], seg2["y1"]
    else:
        x_int2, y_int2 = seg2["x2"], seg2["y2"]

    # Apply axis flipping if necessary
    if args.flip:
        x_int1, y_int1 = y_int1, x_int1
        x_int2, y_int2 = y_int2, x_int2

    # Annotate the internal ends
    plt.text(x_int1, y_int1, f"({int(x_int1)}, {int(y_int1)})", color="blue", fontsize=8, ha="right", va="bottom")
    plt.text(x_int2, y_int2, f"({int(x_int2)}, {int(y_int2)})", color="blue", fontsize=8, ha="right", va="bottom")

base_output = os.path.splitext(args.output)[0]
plot_filename = (base_output + ".pdf") if args.pdf else (base_output + ".png")
plt.savefig(plot_filename)

