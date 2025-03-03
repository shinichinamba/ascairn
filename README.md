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
pip install ascairn (--user)
```

2. Download the resource files for ascairn (in preparation).

# Mini-tutorial

1. Prepare the sequence data

```
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989340/NA12877.final.cram seq_data
```

Alternatively (which might be faster),
```
aws s3 cp s3://1000genomes/1000G_2504_high_coverage/additional_698_related/data/ERR3989340/NA12877.final.cram seq_data/
```

* In fact, if you have access to the AWS S3 BAM file which you want to analyse, and SAMtools is correctly installed, you do not need to download the BAM file and directly specify S3 path.


2. Calculate the sequence depth and determine the sex for the sample.

3. Parse the rare-kmer counts of the target BAM file in the alpha satellite region.

4. 
5. Perform centromere typing on a local BAM/CRAM file:
```
ascairn type /path/to/local/file.bam /path/to/output
```
Example command using a CRAM file located in Amazon S3 (it will take 1~2 hours):
```
ascairn type s3://1000genomes/1000G_2504_high_coverage/data/ERR3240114/HG00096.final.cram output/HG00096 --threads 8
```
