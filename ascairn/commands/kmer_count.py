import click
import os
from ascairn.utils import is_tool, is_exists_bam, convert_tsv_to_fasta, check_kmer_size_from_kmer_fasta, count_rare_kmer

from ascairn.logger import get_logger
logger = get_logger(__name__)

@click.command()
@click.argument("bam_file", type=click.Path(exists=False))
@click.argument("kmer_file", type=click.Path(exists=True))
@click.argument("cen_region_file", type=click.Path(exists=True))
@click.argument("output_file", type=click.Path())
@click.option("-t", "--threads", default=4, help="Number of threads to use.")
def kmer_count_command(bam_file, kmer_file, cen_region_file, output_file, threads):


    # check if the executables exist
    is_tool("samtools")
    is_tool("jellyfish")

    # check input file existences
    is_exists_bam(bam_file)

    # make directory for the output prefix
    output_dir = os.path.dirname(output_file)
    if output_dir != '' and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    is_make_kmer_file_fasta = False
    if kmer_file.endswith(".fa") or kmer_file.endswith(".fasta"):
        logger.info("The provided `kmer_file` has a suffix of '.fa' or '.fasta'; treating it as a FASTA format file.")
        kmer_file_fasta = kmer_file
    else:
        logger.info("Converting the `kmer_file` from TSV format to FASTA format for Jellyfish execution.")
        convert_tsv_to_fasta(kmer_file, output_file + ".tmp.kmer_list.fa")
        kmer_file_fasta = output_file + ".tmp.kmer_list.fa"
        is_make_kmer_file_fasta = True

    kmer_size = check_kmer_size_from_kmer_fasta(kmer_file_fasta)

    count_rare_kmer(bam_file, output_file, cen_region_file, kmer_file_fasta, kmer_size, threads)

    if is_make_kmer_file_fasta:
        os.remove(output_file + ".tmp.kmer_list.fa")

    logger.info("Completed.")


