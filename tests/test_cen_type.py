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



@pytest.mark.parametrize("chrom,is_single_hap", CHROMOSOMES)
def test_cen_type(chrom, is_single_hap, panel_dir, output_dir, expected_dir):
    # Use expected files as input (independent of other tests)
    kmer_count_file = os.path.join(expected_dir, "NA12877.kmer_count.txt")
    depth_file = os.path.join(expected_dir, "NA12877.depth.txt")

    output_prefix = os.path.join(output_dir, f"NA12877.chr{chrom}")

    cmd = [
        "ascairn", "cen_type",
        kmer_count_file,
        "-o", output_prefix,
        "--kmer_info", os.path.join(panel_dir, "kmer_info", f"chr{chrom}.kmer_info.txt.gz"),
        "--hap_info", os.path.join(panel_dir, "hap_info", f"chr{chrom}.hap_info.txt"),
        "--depth_file", depth_file,
    ]
    if is_single_hap:
        cmd.append("--single_hap")

    subprocess.run(cmd, check=True, env=os.environ)

    # Compare cluster and haplotype results with approximate floating point comparison
    for suffix in ["cluster.hap_pair.txt", "haplotype.hap_pair.txt",
                   "cluster.marker_prob.txt", "haplotype.marker_prob.txt"]:
        output_file = f"{output_prefix}.{suffix}"
        expected_file = os.path.join(expected_dir, f"NA12877.chr{chrom}.{suffix}")
        assert os.path.exists(output_file), f"Output file not found: {output_file}"
        assert compare_tsv_approx(output_file, expected_file), f"Mismatch: {suffix}"
