# Test for ascairn somatic_cna command (tiny synthetic data)
import gzip
import subprocess


def test_somatic_cna(tmp_path):
    marker_prob = tmp_path / "marker_prob.txt"
    kmer_info = tmp_path / "chr1.kmer_info.txt.gz"
    normal_count = tmp_path / "normal.kmer_count.txt"
    tumor_count = tmp_path / "tumor.kmer_count.txt"
    out_prefix = tmp_path / "out"
    out = tmp_path / "out.somatic_cna.txt"

    # cols: Marker H1 pos1 cnt1 H2 pos2 cnt2 Prob_00..Prob_22
    # mA -> hap1-specific, mB -> hap2-specific, mC -> ambiguous (dropped),
    # mD -> hap1-specific but no kmer_info pos (dropped), mE -> hap1 with Normal_count 0 (Ratio NA)
    marker_prob.write_text(
        "Marker\tHaplotype1\tMean_marker_pos1\tMarker_count1\tHaplotype2\tMean_marker_pos2\tMarker_count2\t"
        "Prob_00\tProb_01\tProb_02\tProb_10\tProb_11\tProb_12\tProb_20\tProb_21\tProb_22\n"
        "mA\tH1\tNA\tNA\tH2\tNA\tNA\t0.1\t0\t0\t0.9\t0\t0\t0\t0\t0\n"
        "mB\tH1\tNA\tNA\tH2\tNA\tNA\t0.1\t0.9\t0\t0\t0\t0\t0\t0\t0\n"
        "mC\tH1\tNA\tNA\tH2\tNA\tNA\t0.1\t0.45\t0\t0.45\t0\t0\t0\t0\t0\n"
        "mD\tH1\tNA\tNA\tH2\tNA\tNA\t0.1\t0\t0\t0.9\t0\t0\t0\t0\t0\n"
        "mE\tH1\tNA\tNA\tH2\tNA\tNA\t0.1\t0\t0\t0.9\t0\t0\t0\t0\t0\n"
    )

    with gzip.open(kmer_info, "wt") as h:
        h.write(
            "Marker\tHaplotype\tAssembly_source\tContig_len\tMarker_pos\tMarker_strand\n"
            "mA\tH1\tsrc\t1000\t100\t+\n"
            "mB\tH2\tsrc\t1000\t200\t+\n"
            "mE\tH1\tsrc\t1000\t500\t+\n"  # mD intentionally absent on H1
        )

    normal_count.write_text("mA\t10\nmB\t20\nmE\t0\n")
    tumor_count.write_text("mA\t30\nmB\t20\nmE\t5\n")

    subprocess.run(
        ["ascairn", "somatic_cna",
         "--marker_prob", str(marker_prob),
         "--normal_count", str(normal_count),
         "--tumor_count", str(tumor_count),
         "--kmer_info", str(kmer_info),
         "-o", str(out_prefix)],
        check=True,
    )

    lines = out.read_text().splitlines()
    assert lines[0] == "Marker\tHap_pos\tHaplotype\tNormal_count\tTumor_count\tRatio"
    # mC (ambiguous) and mD (no kmer_info pos) are dropped.
    assert lines[1:] == [
        "mA\t100\t1\t10\t30\t3.0",
        "mB\t200\t2\t20\t20\t1.0",
        "mE\t500\t1\t0\t5\tNA",
    ]
