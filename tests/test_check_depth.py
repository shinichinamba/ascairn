# Test for ascairn check_depth command
import os
import subprocess
import filecmp


def test_check_depth(s3_cram, resource_dir, output_dir, expected_dir):
    output_file = os.path.join(output_dir, "NA12877.depth.txt")
    expected_file = os.path.join(expected_dir, "NA12877.depth.txt")

    subprocess.run(
        [
            "ascairn", "check_depth",
            s3_cram,
            os.path.join(resource_dir, "chr22_long_arm_hg38.bed"),
            output_file,
            "--x_region_file", os.path.join(resource_dir, "chrX_short_arm_hg38.bed"),
            "--threads", "4",
        ],
        check=True,
    )

    assert os.path.exists(output_file)
    assert filecmp.cmp(output_file, expected_file, shallow=False)
