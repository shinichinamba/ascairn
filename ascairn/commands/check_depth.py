import click
import os
from ascairn.utils import is_tool, is_exists_bam, check_depth

from ascairn.logger import get_logger
logger = get_logger(__name__)

@click.command()
@click.argument("bam_file", type=click.Path(exists=False))
@click.argument("baseline_region_file", type=click.Path(exists=True))
@click.argument("output_file", type=click.Path())
@click.option("--x_region_file", type=click.Path(exists=True), default = None)
@click.option("--threads", default=4, help="Number of threads to use.")
def check_depth_command(bam_file, baseline_region_file, output_file, x_region_file, threads):

    # check if the executables exist
    is_tool("samtools")
    is_tool("mosdepth")

    # check input file existences
    is_exists_bam(bam_file)

    # make directory for the output prefix
    output_dir = os.path.dirname(output_file)
    if output_dir != '' and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    logger.info("Checking the sequence coverage in the baseline region") 
    depth = check_depth(bam_file, output_file + ".tmp.coverage.txt", baseline_region_file, num_threads = threads)
    os.remove(output_file + ".tmp.coverage.txt") 
    
    if x_region_file is not None:
        logger.info("Checking the sequence coverage in the short arm of chromosome X short-arm region for sex determination")
        depth_chrX = check_depth(bam_file, output_file + ".tmp.chrX_coverage.txt", x_region_file, num_threads = threads)
        os.remove(output_file + ".tmp.chrX_coverage.txt")
    
        logger.info("The depth ratio: %s" % str(depth_chrX / depth))
        if depth_chrX / depth < 0.8: 
            logger.info("Sex has been set to male")
            sex = "male"
        else:
            logger.info("Sex has been set to female")
            sex = "female"

    with open(output_file, 'w') as hout:
        print(f'Coverage: {depth}', file = hout)
        print(f'Baseline region file: {baseline_region_file}', file = hout)
        if x_region_file is not None:
            print(f'Sex: {sex}', file = hout)
            print(f'ChrX region file: {x_region_file}', file = hout)
            print(f'ChrX coverage: {depth_chrX}', file = hout)

    logger.info("Completed.")
    
