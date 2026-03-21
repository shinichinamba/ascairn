import click
import os
from ascairn.match import match_cluster_haplotype, match_cluster_haplotype_single

from ascairn.logger import get_logger
logger = get_logger(__name__)


def read_depth_from_file(depth_file):
    """Parse Coverage value from depth.txt file."""
    with open(depth_file) as f:
        for line in f:
            if line.startswith("Coverage:"):
                return float(line.split()[1])
    raise ValueError(f"Coverage not found in {depth_file}")


@click.command()
@click.argument("kmer_count_file", type=click.Path(exists=True))
@click.option("-o", "--output_prefix", required=True, type=click.Path())
@click.option("--kmer_info", "kmer_info_file", required=True, type=click.Path(exists=True))
@click.option("--hap_info", "hap_info_file", required=True, type=click.Path(exists=True))
@click.option("--depth_file", required=True, type=click.Path(exists=True))
@click.option("--single_hap", is_flag=True, default=False, help="Use single-haplotype mode (e.g. chrX in males).")
@click.option("--hap_candidates", type=click.Path(exists=True), default=None, help="File with haplotype names (one per line) to restrict haplotype-level matching.")
@click.option("--cluster_ratio", default=0.1, hidden=True)
@click.option("--pseudo_count", default=0.1, hidden=True)
@click.option("--nbinom_size_0", default=0.5, hidden=True)
@click.option("--nbinom_size", default=8.0, hidden=True)
@click.option("--nbinom_mu_0_unit", default=0.8 / 30, hidden=True)
@click.option("--nbinom_mu_unit", default=0.4, hidden=True)
def cen_type_command(kmer_count_file, output_prefix, kmer_info_file, hap_info_file, depth_file, single_hap, hap_candidates,
                 cluster_ratio, pseudo_count, nbinom_size_0, nbinom_size, nbinom_mu_0_unit, nbinom_mu_unit):

    depth = read_depth_from_file(depth_file)

    # make directory for the output prefix
    output_dir = os.path.dirname(output_prefix)
    if output_dir != '' and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    logger.info("Matching haplotypes")
    if single_hap:
        match_cluster_haplotype_single(kmer_count_file, output_prefix, kmer_info_file, hap_info_file, depth,
            cluster_ratio=cluster_ratio, pseudo_count=pseudo_count, nbinom_size_0=nbinom_size_0,
            nbinom_size=nbinom_size, nbinom_mu_0_unit=nbinom_mu_0_unit, nbinom_mu_unit=nbinom_mu_unit,
            hap_candidates_file=hap_candidates)
    else:
        match_cluster_haplotype(kmer_count_file, output_prefix, kmer_info_file, hap_info_file, depth,
            cluster_ratio=cluster_ratio, pseudo_count=pseudo_count,
            nbinom_size_0=nbinom_size_0, nbinom_size=nbinom_size, nbinom_mu_0_unit=nbinom_mu_0_unit, nbinom_mu_unit=nbinom_mu_unit,
            hap_candidates_file=hap_candidates)

    logger.info("Completed.")
 
