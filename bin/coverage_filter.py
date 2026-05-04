import gzip
import sys
from collections import Counter

import numpy as np


def weighted_median_mad(x, hist):
    x = np.asarray(x, dtype=np.float64)
    hist = np.asarray(hist, dtype=np.float64)
    order = np.argsort(x)
    x = x[order]
    w = hist[order]
    cdf = np.cumsum(w) / np.sum(w)
    median = np.interp(0.5, cdf, x)
    abs_dev = np.abs(x - median)
    order2 = np.argsort(abs_dev)
    x2 = abs_dev[order2]
    w2 = w[order2]
    cdf2 = np.cumsum(w2) / np.sum(w2)
    mad = np.interp(0.5, cdf2, x2)

    return median, mad


def read_coverage(bed_fn):
    cov_regions = []
    hist_vals = Counter()
    with gzip.open(bed_fn, 'rt') as bed:
        for record in bed:
            contig, start, end, cov = record.split()
            start, end, cov = int(start), int(end), int(cov)
            cov_regions.append((contig, start, end, cov))
            hist_vals[cov] += end - start

    mx = max(hist_vals)
    mn = min(hist_vals)
    hist = np.zeros(mx - mn + 1)
    for i, n in hist_vals.items():
        hist[i - mn] += n
    x = np.log1p(np.arange(mn, mx + 1))
    return cov_regions, x, hist
            

def filter_coverage(cov_regions, min_cov, max_cov):
    region = None
    for contig, start, end, cov in cov_regions:
        if (cov >= min_cov) and (cov <= max_cov):
            if region is None:
                region = [contig, start, end]
            elif contig != region[0]:
                sys.stdout.write(f'{region[0]}\t{region[1]}\t{region[2]}\n')
                region = [contig, start, end]
            else:
                region[2] = end
        else:
            if region is not None:
                sys.stdout.write(f'{region[0]}\t{region[1]}\t{region[2]}\n')
                region = None
    if region is not None:
        sys.stdout.write(f"{region[0]}\t{region[1]}\t{region[2]}\n")


if __name__ == '__main__':
    cov_regions, x, hist = read_coverage(sys.argv[1])
    median, mad = weighted_median_mad(x, hist)
    sigma = 1.4826 * mad
    min_cov = np.expm1(median - 5 * sigma)
    max_cov = np.expm1(median + 5 * sigma)
    sys.stderr.write(f'Min coverage (Median + 5*sigma) = {min_cov}\n')
    sys.stderr.write(f'Max coverage (Median + 5*sigma) = {max_cov}\n')
    filter_coverage(cov_regions, min_cov=min_cov, max_cov=max_cov)