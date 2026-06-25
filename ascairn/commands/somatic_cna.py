import os
import gzip
import csv
import click

# Build the D_M-equivalent somatic CNA marker table for one centromere.
# For each marker that ascairn cen_type assigned to one proxy haplotype (hap1 or hap2),
# look up its position on that haplotype from the panel kmer_info and pair the
# normal/tumor rare k-mer counts. No imputation, no filtering, no binning.


def _open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def _read_best_haps(marker_prob_file):
    # The best hap pair is encoded in marker_prob.txt's Haplotype1/Haplotype2 columns,
    # but markers absent from a hap have "NA" there, and the two columns are filled
    # independently. Scan each column for its first non-NA value.
    hap1 = hap2 = None
    with open(marker_prob_file) as hin:
        for F in csv.DictReader(hin, delimiter="\t"):
            if hap1 is None and F["Haplotype1"] != "NA":
                hap1 = F["Haplotype1"]
            if hap2 is None and F["Haplotype2"] != "NA":
                hap2 = F["Haplotype2"]
            if hap1 is not None and hap2 is not None:
                break
    if hap1 is None or hap2 is None:
        raise ValueError(f"Could not determine best haplotypes from {marker_prob_file}")
    return hap1, hap2


def _build_hap_pos(kmer_info_file, hap1, hap2):
    # One pass over the (large) kmer_info, collecting absolute Marker_pos per
    # (marker, hap) for the two proxy haplotypes only. Panel Marker_pos already has
    # the 100kb flank encoded as negative values, so Hap_pos = Marker_pos directly
    # (no -100000 offset; verified vs Mean_marker_pos1).
    # Contig_len is constant within a haplotype; grab it for the plot bin width.
    marker_hap_pos_raw = {}
    hap_len = {}
    with _open_maybe_gzip(kmer_info_file) as hin:
        for F in csv.DictReader(hin, delimiter="\t"):
            hap = F["Haplotype"]
            if hap != hap1 and hap != hap2:
                continue
            marker_hap_pos_raw.setdefault((F["Marker"], hap), []).append(int(F["Marker_pos"]))
            if hap not in hap_len:
                hap_len[hap] = int(F["Contig_len"])

    # Aggregate a single haplotype position per (marker, hap). Same rule as
    # hap_pos_impute.py: one occurrence -> use it; two within 100kb -> average;
    # otherwise ambiguous -> drop.
    marker_hap_pos = {}
    for key, pos_vec in marker_hap_pos_raw.items():
        if len(pos_vec) == 1:
            marker_hap_pos[key] = pos_vec[0]
        elif len(pos_vec) == 2:
            if abs(pos_vec[0] - pos_vec[1]) < 100000:
                marker_hap_pos[key] = (pos_vec[0] + pos_vec[1]) / 2
        # len > 2: leave unassigned (ambiguous)

    return marker_hap_pos, hap_len


def _read_counts(count_file):
    # kmer_count.txt is headerless: Marker<TAB>count
    counts = {}
    with open(count_file) as hin:
        for line in hin:
            marker, count = line.rstrip("\n").split("\t")
            counts[marker] = int(count)
    return counts


def _read_clusters(hap_info_file, hap1, hap2):
    # hap_info maps Haplotype -> Cluster (the aHOR-HG). Pull it for the two best haps.
    cluster = {}
    with open(hap_info_file) as hin:
        for F in csv.DictReader(hin, delimiter="\t"):
            if F["Haplotype"] == hap1:
                cluster[1] = F["Cluster"]
            elif F["Haplotype"] == hap2:
                cluster[2] = F["Cluster"]
    return cluster


def _plot_somatic_cna(records, hap_len, hap_name, hap_cluster, output_pdf, title):
    # Per-bin Tumor/Normal ratio boxplots, one panel per proxy haplotype.
    # Mirrors plot_somatic_cna.R: 20 bins over the true haplotype length, same QC.
    # records: list of (hap_num, hap_pos, normal_count, ratio). Raises ImportError
    # if matplotlib is missing (caller warns and skips).
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    def bin_edges(length):
        bin_size = (length / 1e6) / 20
        return np.arange(-bin_size, (length / 1e6) + 2 * bin_size, bin_size)

    def panel_data(hap_num):
        length = hap_len.get(hap_num)
        if length is None:
            return None
        edges = bin_edges(length)
        # QC: require Normal_count >= 8 (drop low-depth markers).
        xs, ys = [], []
        for hn, pos, normal, ratio in records:
            if hn != hap_num or ratio == "NA":
                continue
            if normal < 8:
                continue
            xs.append(pos / 1e6)
            ys.append(ratio)
        if not ys:
            return None
        # idx in [0, len(edges)-2] selects the bin; drop points outside the range.
        bins = [[] for _ in range(len(edges) - 1)]
        idx = np.digitize(xs, edges) - 1
        for k, b in enumerate(idx):
            if 0 <= b < len(bins):
                bins[b].append(ys[k])
        return edges, bins

    panels = [(h, panel_data(h)) for h in (1, 2)]
    panels = [(h, d) for h, d in panels if d is not None]
    if not panels:
        return False

    fig, axes = plt.subplots(len(panels), 1, figsize=(12 / 2.54, 5 / 2.54 * len(panels)),
                             squeeze=False)
    for ax, (hap_num, (edges, bins)) in zip(axes[:, 0], panels):
        positions = [i + 1 for i in range(len(bins))]
        data = [b if b else [np.nan] for b in bins]
        bp = ax.boxplot(data, positions=positions, widths=0.6, patch_artist=True,
                        showfliers=True, flierprops=dict(markersize=1),
                        whiskerprops=dict(linewidth=0.2), capprops=dict(linewidth=0.2))
        for box in bp["boxes"]:
            box.set(facecolor="#FFBC42", linewidth=0.4)
        # R cut()-style interval labels "(lo,hi]" with fixed 2 decimals.
        labels = [f"({edges[i]:.2f},{edges[i + 1]:.2f}]" for i in range(len(edges) - 1)]
        ax.set_xticks(positions)
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=5,
                           rotation_mode="anchor")
        ax.set_xlabel("aHOR-hap position bin (Mb)", fontsize=6)
        ax.set_ylabel("Rare k-mer count ratio", fontsize=6)
        sub = f"Proxy aHOR-hap {hap_num}"
        nm = hap_name.get(hap_num)
        if nm:
            sub += f": {nm}"
        cl = hap_cluster.get(hap_num)
        if cl is not None:
            sub += f" (aHOR-HG {cl})"
        ax.set_title(f"{title}; {sub}", fontsize=6)
        ax.tick_params(labelsize=5)
        # Pull the rotated x tick labels close to the axis (mirrors R's margin t=-0.5).
        ax.tick_params(axis="x", pad=0)

    fig.tight_layout()
    fig.savefig(output_pdf)
    plt.close(fig)
    return True


@click.command()
@click.option("--marker_prob", "marker_prob_file", required=True, type=click.Path(exists=True),
              help="Normal cen_type *.haplotype.marker_prob.txt")
@click.option("--normal_count", "normal_count_file", required=True, type=click.Path(exists=True),
              help="Normal kmer_count.txt (headerless: Marker, count)")
@click.option("--tumor_count", "tumor_count_file", required=True, type=click.Path(exists=True),
              help="Tumor kmer_count.txt (headerless: Marker, count)")
@click.option("--kmer_info", "kmer_info_file", required=True, type=click.Path(exists=True),
              help="Panel chr{N}.kmer_info.txt.gz")
@click.option("--hap_info", "hap_info_file", default=None, type=click.Path(exists=True),
              help="Panel chr{N}.hap_info.txt; adds the aHOR-HG (cluster) to plot titles.")
@click.option("-o", "--output_prefix", required=True, type=click.Path())
@click.option("--prob_hap_thres", default=0.8, hidden=True)
@click.option("--prob_other_thres", default=0.05, hidden=True)
def somatic_cna_command(marker_prob_file, normal_count_file, tumor_count_file,
                        kmer_info_file, hap_info_file, output_prefix,
                        prob_hap_thres, prob_other_thres):

    output_file = output_prefix + ".somatic_cna.txt"

    hap1, hap2 = _read_best_haps(marker_prob_file)
    marker_hap_pos, hap_len_by_name = _build_hap_pos(kmer_info_file, hap1, hap2)

    normal_counts = _read_counts(normal_count_file)
    tumor_counts = _read_counts(tumor_count_file)

    # Collected for the plot: (hap_num, hap_pos, normal_count, ratio).
    records = []

    with open(output_file, "w") as hout:
        print("Marker\tHap_pos\tHaplotype\tNormal_count\tTumor_count\tRatio", file=hout)

        with open(marker_prob_file) as hin:
            for F in csv.DictReader(hin, delimiter="\t"):
                marker = F["Marker"]
                p_hap1 = float(F["Prob_10"]) + float(F["Prob_20"])
                p_hap2 = float(F["Prob_01"]) + float(F["Prob_02"])

                if p_hap1 >= prob_hap_thres and p_hap2 <= prob_other_thres:
                    hap_num, hap = 1, hap1
                elif p_hap2 >= prob_hap_thres and p_hap1 <= prob_other_thres:
                    hap_num, hap = 2, hap2
                else:
                    continue

                # Observed markers only (no imputation).
                if (marker, hap) not in marker_hap_pos:
                    continue

                hap_pos = marker_hap_pos[(marker, hap)]
                normal_count = normal_counts.get(marker, 0)
                tumor_count = tumor_counts.get(marker, 0)
                ratio = tumor_count / normal_count if normal_count > 0 else "NA"

                print(f"{marker}\t{hap_pos}\t{hap_num}\t{normal_count}\t{tumor_count}\t{ratio}", file=hout)
                records.append((hap_num, hap_pos, normal_count, ratio))

    # Always attempt the plot; skip with a warning if matplotlib is unavailable.
    output_pdf = output_prefix + ".somatic_cna.pdf"
    hap_len = {1: hap_len_by_name.get(hap1), 2: hap_len_by_name.get(hap2)}
    hap_name = {1: hap1, 2: hap2}
    hap_cluster = _read_clusters(hap_info_file, hap1, hap2) if hap_info_file else {}
    title = os.path.basename(output_prefix)
    try:
        drawn = _plot_somatic_cna(records, hap_len, hap_name, hap_cluster, output_pdf, title)
        if not drawn:
            click.echo(f"No data passed QC; skipped plot {output_pdf}", err=True)
    except ImportError:
        click.echo("matplotlib not installed; skipping plot "
                    "(install with: pip install ascairn[plot])", err=True)
