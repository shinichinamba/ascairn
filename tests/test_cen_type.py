# Test for ascairn cen_type command
import os
import subprocess
import pytest


def compare_tsv_approx(file1, file2, rel_tol=1e-9):
    """Compare two TSV files with approximate floating point comparison."""
    with open(file1) as f1, open(file2) as f2:
        for line1, line2 in zip(f1, f2):
            cols1 = line1.strip().split('\t')
            cols2 = line2.strip().split('\t')
            if len(cols1) != len(cols2):
                return False
            for c1, c2 in zip(cols1, cols2):
                if c1 == c2:
                    continue
                try:
                    v1, v2 = float(c1), float(c2)
                    if abs(v1 - v2) > rel_tol * max(abs(v1), abs(v2), 1e-10):
                        return False
                except ValueError:
                    return False
    return True

# Test parameters: (chromosome, is_single_hap)
# NA12877 is male, so chrX uses --is_single_hap
# Note: chr1 excluded due to high memory requirement (~6GB)
CHROMOSOMES = [
    ("22", False),
    ("X", True),
]


def read_depth_from_file(depth_file):
    with open(depth_file) as f:
        for line in f:
            if line.startswith("Coverage:"):
                return line.split()[1]
    raise ValueError("Coverage not found in depth file")


@pytest.mark.parametrize("chrom,is_single_hap", CHROMOSOMES)
def test_cen_type(chrom, is_single_hap, resource_dir, resource_version, output_dir, expected_dir):
    # Use expected files as input (independent of other tests)
    kmer_count_file = os.path.join(expected_dir, "NA12877.kmer_count.txt")
    depth_file = os.path.join(expected_dir, "NA12877.depth.txt")
    depth = read_depth_from_file(depth_file)

    output_prefix = os.path.join(output_dir, f"NA12877.chr{chrom}")

    # cluster directory name differs by version
    cluster_dir = "cluster_m3" if resource_version == "ver_2024-12-06" else "cluster"

    cmd = [
        "ascairn", "cen_type",
        kmer_count_file,
        output_prefix,
        os.path.join(resource_dir, "kmer_info", f"chr{chrom}.kmer_info.txt.gz"),
        os.path.join(resource_dir, cluster_dir, f"chr{chrom}.cluster_marker_count.txt.gz"),
        depth,
        "--cluster_haplotype_file", os.path.join(resource_dir, cluster_dir, f"chr{chrom}.hap_cluster.txt"),
    ]
    if is_single_hap:
        cmd.append("--is_single_hap")

    subprocess.run(cmd, check=True, env=os.environ)

    # Compare cluster and haplotype results with approximate floating point comparison
    for suffix in ["cluster.hap_pair.txt", "haplotype.hap_pair.txt",
                   "cluster.marker_prob.txt", "haplotype.marker_prob.txt"]:
        output_file = f"{output_prefix}.{suffix}"
        expected_file = os.path.join(expected_dir, f"NA12877.chr{chrom}.{suffix}")
        assert os.path.exists(output_file), f"Output file not found: {output_file}"
        assert compare_tsv_approx(output_file, expected_file), f"Mismatch: {suffix}"
