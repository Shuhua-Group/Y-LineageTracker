import time
import logging
import argparse

from ClassifyHaplogroup import MatchVariant


'''
Description:
This module is used for genotyping Y-STRs from BAM or VCF indels.
Different Y-STR panels are supported in this command and can be selected by the –-panel option.
The accuracy and genotyping rate depend on the quality of input data (such as reads coverage of BAM files).
'''


def genostr_parser():

    parser = argparse.ArgumentParser('genostr', description='(c) Y-LineageTracker: Y-STR genotyping')
    # function used for Y-STR genotyping
    parser.add_argument('genostr',
                        help='Genotype Y-STR from BAM or VCF indels.')
    # required, input file, BAM or VCF
    input = parser.add_mutually_exclusive_group(required=True)
    input.add_argument('--vcf',
                        type=str,
                        action='store',
                        help='vcf: Take VCF as input file, only one VCF is required.')
    input.add_argument('--bam',
                        type=str,
                        action='append',
                        help='bam: Take BAMs or a list of BAM files for analysis')
    input.add_argument('--cram',
                        type=str,
                        action='append',
                        help='cram: Take CRAMs or a list of BAM files for analysis')
    # optional, provide samples to be analyzed
    parser.add_argument('-s', '--sample',
                        required=False,
                        type=str,
                        action='store',
                        help='Sample: A file containing sample individuals only.')
    # optional, version of reference genome: 37 or 38
    parser.add_argument('-b', '--build',
                        required=True,
                        type=int,
                        action='store',
                        choices=[37, 38],
                        help='build: The build version of reference genome used for Y-STR genotyping. (37 or 38).')
    # optional, filter females by misisng rate
    parser.add_argument('--filter',
                        required=False,
                        type=str,
                        action='store',
                        nargs='?',
                        const='0.75',
                        help='filter: Filter female individuals by setting cutoff of missing rate, default is 0.75. If the missing rate value is not set, the program will automatically filter female individuals by the built-in algorithm.')
    # optional, the set of Y-STR loci used for genotyping
    parser.add_argument('--panel',
                        required=False,
                        type=str,
                        default='named',
                        choices=['all', 'named', 'minimal', 'ppy', 'pp23', 'yf'],
                        help='panel: The set of Y-STR loci for genotyping. These sets are referenced from commonly used STR multiplex kits:\n'
                             'all: all the possible Y-STR loci\n'
                             'named: all the possible named Y-STR loci\n'
                             'minimal: minimal panel for Y-STR loci\n'
                             'ppy: PowerPlex Y panel for Y-STR loci\n'
                             'pp23: PowerPlex 23 panel for Y-STR loci\n'
                             'Yf: Yfiler panel for Y-STR loci')
    # optional, the prefix of output
    parser.add_argument('-o', '--output',
                        required=False,
                        type=str,
                        action='store',
                        help='output: The prefix of output files.')

    args = parser.parse_args()

    return args


# print program information and write to log file
def set_log(input, log_file, format, cutoff, args_log):

    logger = logging.getLogger()
    logger.setLevel(level=logging.INFO)

    handler = logging.FileHandler(log_file, mode='w')
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter('[%(asctime)s] - [%(levelname)s]: %(message)s')
    handler.setFormatter(formatter)

    console = logging.StreamHandler()
    console.setLevel(logging.INFO)

    logger.addHandler(handler)
    logger.addHandler(console)

    log_info = ['[Y-LineageTracker] [GenotypeSTR]',
                '[Y-LineageTracker] Run Date: ' + time.asctime(time.localtime(time.time())),
                '[Y-LineageTracker] Input File: %s' % input,
                '[Y-LineageTracker] File Format: %s' % format,
                '[Y-LineageTracker] Reference Build: GRCh%d' % args_log.build,
                '[Y-LineageTracker] STR Panel: %s' % args_log.panel]

    if args_log.sample:
        log_info.append('[Y-LineageTracker] List of samples provided')
    else:
        if not cutoff is None:
            if cutoff:
                log_info.append('[Y-LineageTracker] Missing rate cutoff: %.2f' % cutoff)
            else:
                log_info.append('[Y-LineageTracker] Filter female individuals automatically')

    print('\n')
    for i in log_info:
        logger.info(i)

# names of output files
def check_str_output(output):

    from FilesIO import get_out_path
    path = get_out_path(output, 'Genostr')
    out_file = path + '.str.matrix'
    log_file = path + '.StrGenotype.log'

    return out_file, log_file


def main():

    start = time.perf_counter()
    arguments = genostr_parser()

    from FilesIO import check_cutoff, check_BAM_VCF_input, check_samples, time_count
    # check input files
    input, format = check_BAM_VCF_input(arguments.vcf, arguments.bam, arguments.cram)

    # check cutoff: cutoff is a float between 0 and 1, else False
    cutoff = check_cutoff(arguments.filter, format)

    # set of output files
    out_file, log_file = check_str_output(arguments.output)

    # set log file
    set_log(input, log_file, format, cutoff, arguments)

    # check whether provide sample list
    samples = check_samples(arguments.sample)

    # filter samples
    from ProcessData import ReadFilterData
    read_filter = ReadFilterData.ReadFilterData(input, samples)

    # if BAM, genotype Y-STR from ReadFilterData module
    if format == 'bam':
        result = read_filter.read_bam(arguments.panel, arguments.build, 'str', format)[0]
    elif format == 'cram':
        result = read_filter.read_bam(arguments.panel, arguments.build, 'str', format)[0]
    # if VCF filter samples from ReadFilterData, then genotype Y-STRs from MatchVariant
    else:
        data_info, inds = read_filter.read_filter_variant(cutoff, format)
        # genotype Y-STR
        match = MatchVariant.MatchVariant(arguments.build)
        result = match.match_variant_str(data_info, inds, arguments.panel)

    # output result in matrix format
    result.index.name = 'SampleID'
    result.to_csv(out_file,
                  sep='\t',
                  index=True)

    time_count(start)


if __name__ == '__main__':
    main()
