<div align="center">
  <img src="image/karamatsu_cairn.png" alt="Cairn at Mt. Karamatsu" width="250">
  <p>Maruyama cairn at Mt. Karamatsu</p>
</div>

# ascairn

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![CI](https://github.com/friend1ws/ascairn/actions/workflows/python-test.yml/badge.svg)

`ascairn` (alpha-satellite cairn) is software for estimating centromere variation from short-read sequencing data using rare k-mers within centromere sequences.
For each chromosome, ascairn identifies the most likely pair of centromere clusters and the nearest haplotypes from a reference panel ([Shiraishi et al., bioRxiv, 2025](https://doi.org/10.1101/2025.07.26.666712)).

<div align="center">
  <img src="image/categorization_overview.png" alt="Overview of centromere classification using rare k-mer" width="750">
</div>

## Prerequisites

### Software
- [samtools](https://github.com/samtools/samtools)
- [Jellyfish](https://github.com/gmarcais/Jellyfish)
- [mosdepth](https://github.com/brentp/mosdepth)

### Python packages
- click
- scipy
- polars
- boto3 (required only for accessing CRAM files in Amazon S3)

## Installation

1. Install prerequisite software and ensure they are accessible via your `PATH`.

2. Install `ascairn`.
```bash
git clone https://github.com/friend1ws/ascairn.git
cd ascairn
pip install .   # or pip install -e . for development
```

3. Download [ascairn resource files](https://github.com/friend1ws/ascairn_resource).
This repository contains the reference data required by ascairn, including:
   - Rare k-mer list for alpha satellite sequences (`rare_kmer_list.fa`)
   - Centromere region BED files for hg38 and chm13
   - Per-chromosome k-mer information (`kmer_info/`)
   - Per-chromosome haplotype-to-cluster mapping (`hap_info/`)

```bash
git clone https://github.com/friend1ws/ascairn_resource.git
```

After installation, your directory structure should look like this:
```
ascairn/
ascairn_resource/
└── resource/
    ├── common/
    │   ├── cen_region_hg38.bed
    │   ├── cen_region_chm13.bed
    │   ├── chr22_long_arm_hg38.bed
    │   ├── chr22_long_arm_chm13.bed
    │   ├── chrX_short_arm_hg38.bed
    │   └── chrX_short_arm_chm13.bed
    └── panel/
        └── ascairn_paper_2025/
            ├── rare_kmer_list.fa
            ├── kmer_info/
            │   ├── chr1.kmer_info.txt.gz
            │   ├── ...
            │   └── chrX.kmer_info.txt.gz
            └── hap_info/
                ├── chr1.hap_info.txt
                ├── ...
                └── chrX.hap_info.txt
```

## Quick Start

The simplest way to run ascairn is through the wrapper script `ascairn_type_allchr.sh`,
which executes all three steps (`check_depth`, `kmer_count`, `cen_type`) sequentially for chr1-22 and chrX.

### 1. Prepare the sequence data

Download sequence data (aligned to the GRCh38 reference genome) using one of the following methods:

**Option 1: AWS S3 (recommended for faster downloads)**
```
aws s3 cp s3://1000genomes/1000G_2504_high_coverage/additional_698_related/data/ERR3989340/NA12877.final.cram seq_data/
aws s3 cp s3://1000genomes/1000G_2504_high_coverage/additional_698_related/data/ERR3989340/NA12877.final.cram.crai seq_data/
```

> If you have direct access to AWS S3 BAM or CRAM files and SAMtools is properly installed, you can skip downloading the files locally. Instead, directly specify the S3 path in subsequent steps.

**Option 2: FTP**
```
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989340/NA12877.final.cram -P seq_data/
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989340/NA12877.final.cram.crai -P seq_data/
```

### 2. Run the workflow

```bash
bash ascairn_type_allchr.sh seq_data/NA12877.final.cram output/NA12877 ascairn_resource/resource/panel/ascairn_paper_2025 hg38 8
```

**Usage:**
```
ascairn_type_allchr.sh <BAM_FILE> <OUTPUT_PREFIX> <DATA_DIR> <REFERENCE: hg38|chm13> [THREAD_NUM]
```

| Argument         | Description                                              | Default |
|------------------|----------------------------------------------------------|---------|
| First argument   | Path to BAM or CRAM file                                 | -       |
| Second argument  | Output path prefix                                       | -       |
| Third argument   | Path to ascairn resource panel directory                  | -       |
| Fourth argument  | Reference genome for BAM file (`hg38` or `chm13`)       | -       |
| Fifth argument   | Number of threads to use                                 | 8       |

### 3. Check the results

After successful execution, the main output file is:
```
output/NA12877.cen_type.result.txt
```

See [Output Format](#output-format) for how to interpret the results.

## Commands

The workflow consists of three commands that are normally run in order.
You can also run them individually for more control.

### `check_depth`

Checks sequence coverage in a reference region (chr22 long arm) and determines biological sex by comparing chrX coverage.

```bash
ascairn check_depth \
    seq_data/NA12877.final.cram \
    -o output/NA12877.depth.txt \
    --baseline_region ascairn_resource/resource/common/chr22_long_arm_hg38.bed \
    --x_region ascairn_resource/resource/common/chrX_short_arm_hg38.bed \
    -t 8
```

**Output** (`output/NA12877.depth.txt`):
```
Coverage: 37.24
Baseline region file: ascairn_resource/resource/ver_2024-12-06/chr22_long_arm_hg38.bed
Sex: male
ChrX region file: ascairn_resource/resource/ver_2024-12-06/chrX_short_arm_hg38.bed
ChrX coverage: 18.0
```

### `kmer_count`

Extracts reads aligned to alpha satellite regions and counts occurrences of predefined rare k-mers.

```bash
ascairn kmer_count \
    seq_data/NA12877.final.cram \
    -o output/NA12877.kmer_count.txt \
    --kmer_file ascairn_resource/resource/panel/ascairn_paper_2025/rare_kmer_list.fa \
    --cen_region ascairn_resource/resource/common/cen_region_hg38.bed \
    -t 8
```

**Output** (`output/NA12877.kmer_count.txt`): a two-column TSV with k-mer sequences and their counts.

### `cen_type`

Identifies the most likely centromere cluster pair and nearest haplotype pair for a given chromosome. Requires the depth file from `check_depth` output.

```bash
ascairn cen_type \
    output/NA12877.kmer_count.txt \
    -o output/NA12877.chr22 \
    --kmer_info ascairn_resource/resource/panel/ascairn_paper_2025/kmer_info/chr22.kmer_info.txt.gz \
    --hap_info ascairn_resource/resource/panel/ascairn_paper_2025/hap_info/chr22.hap_info.txt \
    --depth_file output/NA12877.depth.txt
```

For male samples on chrX, add the `--single_hap` option (the wrapper script handles this automatically based on `check_depth` output).

**Output files** (using `output/NA12877.chr22` as the output prefix):

| File | Description |
|------|-------------|
| `*.cluster.hap_pair.txt` | Ranked cluster pairs with log-likelihoods. The first data row is the best pair. |
| `*.haplotype.hap_pair.txt` | Ranked haplotype pairs with log-likelihoods. The first data row is the best pair. |
| `*.cluster.marker_prob.txt` | Per-marker copy number probabilities used for cluster assignment |
| `*.haplotype.marker_prob.txt` | Per-marker copy number probabilities used for haplotype assignment |

## Output Format

### `cen_type.result.txt`

The main result file (generated by the wrapper script) is a TSV with one row per chromosome:

| Column | Description |
|--------|-------------|
| Chr | Chromosome (chr1-chr22, chrX) |
| Cluster_1 | Best-matching cluster ID for the first allele |
| Cluster_2 | Best-matching cluster ID for the second allele |
| Haplotype_1 | Nearest haplotype ID for the first allele (e.g., `HG00268.mat`) |
| Haplotype_2 | Nearest haplotype ID for the second allele |

For male chrX, Cluster_2 and Haplotype_2 are `NA`.

### Cluster and haplotype pair files

These files list all candidate pairs ranked by log-likelihood (higher = better fit):

| Column | Description |
|--------|-------------|
| Cluster1 / Haplotype1 | First member of the pair |
| Cluster2 / Haplotype2 | Second member of the pair |
| Loglikelihood | Log-likelihood of the pair given the observed k-mer counts |

## Citation

Shiraishi Y, Ochi Y, Sugawa M, Sakamoto Y, Kimura K, Tsujimura T, Okada A, Okuda R, Namba S, Miyauchi T, Mateos RN, Suzuki H, Chiba K, Ito Y, Nakamura W, Ohka F, Motomura K, Yamamoto T, Kawai Y, Okada Y, Suzuki H, Kato M, Saito R, Garrison E, Logsdon GA, Ogawa S. Rare k-mers reveal centromere haplogroups underlying human diversity and cancer translocations. *bioRxiv*. 2025. doi: [10.1101/2025.07.26.666712](https://doi.org/10.1101/2025.07.26.666712)
