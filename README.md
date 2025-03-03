<div align="center">
  <img src="image/karamatsu_cairn.png" alt="Cairn at Mt. Karamatsu" width="250">
  <p>Maruyama cairn at Mt. Karamatsu</p>
</div>

# ascairn
Centromere sequence analysis using rare k-mer markers

# Dependency
## Software
- [samtools](https://github.com/samtools/samtools)
- [Jellyfish](https://github.com/gmarcais/Jellyfish)
- [mosdepth](https://github.com/brentp/mosdepth)

## Python
- click
- boto3 (required only for accessing CRAM files in Amazon S3)
  
# Install

1. Install all the prerequisite software and install `ascairn`.
```
git clone https://github.com/friend1ws/ascairn.git
cd ascairn
pip install . (--user)
```

2. Download the resource files for ascairn.
```
git clone https://github.com/friend1ws/ascairn_data.git
```
# Mini-tutorial

Using the `ascairn_type_allchr.sh` script, you can execute a standard ascairn workflow.

## Workflow 
- `check_depth`　
  - Check the sequence coverage focusing on a reference region (long arm of chromosome 22)
  - Determine biological sex by assessing the coverage within a specified region of chromosome X (restricted here to the short arm) and take the ratio to the reference region. 

- `kmer_count` 
  - Extract reads aligned to alpha satellite regions and count the number of pre-defined rare k-mers.

- `type`
  - Identify centromeric cluster pairs and the closest haplotype pairs for chromosomes 1–22 and X.

## Step

### 1. Prepare the sequence data

Download the sequence data using one of the following methods:

**Option 1: AWS S3 (recommended for faster downloads)**
```
aws s3 cp s3://1000genomes/1000G_2504_high_coverage/additional_698_related/data/ERR3989340/NA12877.final.cram seq_data/
aws s3 cp s3://1000genomes/1000G_2504_high_coverage/additional_698_related/data/ERR3989340/NA12877.final.cram.crai seq_data/
```

>  If you have direct access to AWS S3 BAM files and SAMtools is properly installed, you can skip downloading the file locally. Instead, directly specify the S3 path to your BAM or CRAM file in subsequent steps.


**Option 2: FTP**
```
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989340/NA12877.final.cram seq_data
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989340/NA12877.final.cram.crai seq_data
```


### 2. Execute the `ascairn_type_allchr.sh` script
Run the following command (runtime: approximately 20–30 minutes):
```
bash ascairn_type_allchr.sh seq_data/NA12877.final.cram output/NA12877 ascairn_data 8
```

**Argument Descriptions:**
- First argument: Path to BAM or CRAM file.
- Second argument: Path to the output prefix.
- Third argument: Path to the ascairn resource data.
- Fourth argument: Number of threads to use (optional, default: 8).

**Result**
After successful execution, you will find the output file at:
```
output/NA12877.cen_type.result.txt
```

