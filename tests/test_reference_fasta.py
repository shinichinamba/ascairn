# Test that -r/--reference_fasta lets ascairn process a CRAM when
# REF_PATH/REF_CACHE are not available for reference resolution.
import os
import subprocess

from test_check_depth import compare_depth_values


def _env_without_ref():
    env = os.environ.copy()
    env.pop("REF_PATH", None)
    env.pop("REF_CACHE", None)
    return env


def test_check_depth_with_reference_fasta(
    s3_cram, reference_fasta, common_dir, output_dir, expected_dir
):
    output_file = os.path.join(output_dir, "NA12877.with_ref_fa.depth.txt")
    expected_file = os.path.join(expected_dir, "NA12877.depth.txt")

    subprocess.run(
        [
            "ascairn", "check_depth",
            s3_cram,
            "-o", output_file,
            "--baseline_region", os.path.join(common_dir, "chr22_long_arm_hg38.bed"),
            "--x_region", os.path.join(common_dir, "chrX_short_arm_hg38.bed"),
            "-r", reference_fasta,
            "-t", "4",
        ],
        env=_env_without_ref(),
        check=True,
    )

    assert os.path.exists(output_file)
    assert compare_depth_values(output_file, expected_file)


def test_kmer_count_with_reference_fasta(
    s3_cram, reference_fasta, common_dir, panel_dir, output_dir, expected_dir
):
    output_file = os.path.join(output_dir, "NA12877.with_ref_fa.kmer_count.txt")
    expected_file = os.path.join(expected_dir, "NA12877.kmer_count.txt")

    subprocess.run(
        [
            "ascairn", "kmer_count",
            s3_cram,
            "-o", output_file,
            "--kmer_file", os.path.join(panel_dir, "rare_kmer_list.fa"),
            "--cen_region", os.path.join(common_dir, "cen_region_curated_margin_hg38.bed"),
            "-r", reference_fasta,
            "-t", "4",
        ],
        env=_env_without_ref(),
        check=True,
    )

    assert os.path.exists(output_file)
    with open(output_file) as f1, open(expected_file) as f2:
        assert f1.read() == f2.read()
