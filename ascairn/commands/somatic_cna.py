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
    marker_hap_pos_raw = {}
    with _open_maybe_gzip(kmer_info_file) as hin:
        for F in csv.DictReader(hin, delimiter="\t"):
            hap = F["Haplotype"]
            if hap != hap1 and hap != hap2:
                continue
            marker_hap_pos_raw.setdefault((F["Marker"], hap), []).append(int(F["Marker_pos"]))

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

    return marker_hap_pos


def _read_counts(count_file):
    # kmer_count.txt is headerless: Marker<TAB>count
    counts = {}
    with open(count_file) as hin:
        for line in hin:
            marker, count = line.rstrip("\n").split("\t")
            counts[marker] = int(count)
    return counts


@click.command()
@click.option("--marker_prob", "marker_prob_file", required=True, type=click.Path(exists=True),
              help="Normal cen_type *.haplotype.marker_prob.txt")
@click.option("--normal_count", "normal_count_file", required=True, type=click.Path(exists=True),
              help="Normal kmer_count.txt (headerless: Marker, count)")
@click.option("--tumor_count", "tumor_count_file", required=True, type=click.Path(exists=True),
              help="Tumor kmer_count.txt (headerless: Marker, count)")
@click.option("--kmer_info", "kmer_info_file", required=True, type=click.Path(exists=True),
              help="Panel chr{N}.kmer_info.txt.gz")
@click.option("-o", "--output_file", required=True, type=click.Path())
@click.option("--prob_hap_thres", default=0.8, hidden=True)
@click.option("--prob_other_thres", default=0.05, hidden=True)
def somatic_cna_command(marker_prob_file, normal_count_file, tumor_count_file,
                        kmer_info_file, output_file,
                        prob_hap_thres, prob_other_thres):

    hap1, hap2 = _read_best_haps(marker_prob_file)
    marker_hap_pos = _build_hap_pos(kmer_info_file, hap1, hap2)

    normal_counts = _read_counts(normal_count_file)
    tumor_counts = _read_counts(tumor_count_file)

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
