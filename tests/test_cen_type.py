# Test for ascairn cen_type command
import os
import subprocess
import filecmp
import pytest

# Test parameters: (chromosome, is_single_hap)
# NA12877 is male, so chrX uses --is_single_hap
CHROMOSOMES = [
    ("1", False),
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
def test_cen_type(chrom, is_single_hap, resource_dir, output_dir, expected_dir):
    # Use expected files as input (independent of other tests)
    kmer_count_file = os.path.join(expected_dir, "NA12877.kmer_count.txt")
    depth_file = os.path.join(expected_dir, "NA12877.depth.txt")
    depth = read_depth_from_file(depth_file)

    output_prefix = os.path.join(output_dir, f"NA12877.chr{chrom}")

    cmd = [
        "ascairn", "cen_type",
        kmer_count_file,
        output_prefix,
        os.path.join(resource_dir, "kmer_info", f"chr{chrom}.kmer_info.txt.gz"),
        os.path.join(resource_dir, "cluster_m3", f"chr{chrom}.cluster_marker_count.txt.gz"),
        depth,
        "--cluster_haplotype_file", os.path.join(resource_dir, "cluster_m3", f"chr{chrom}.hap_cluster.txt"),
    ]
    if is_single_hap:
        cmd.append("--is_single_hap")

    subprocess.run(cmd, check=True)

    # Compare cluster and haplotype results
    for suffix in ["cluster.hap_pair.txt", "haplotype.hap_pair.txt"]:
        output_file = f"{output_prefix}.{suffix}"
        expected_file = os.path.join(expected_dir, f"NA12877.chr{chrom}.{suffix}")
        assert os.path.exists(output_file)
        assert filecmp.cmp(output_file, expected_file, shallow=False)
