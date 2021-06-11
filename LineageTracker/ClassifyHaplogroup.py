import time
import logging
import argparse

from GetLineage import MatchVariant, TrackLineage


'''
Description:
This module is used for NRY haplogroup classification.
The program will classify NRY haplogroup and construct the lineage of each individual based on the ISOGG NRY tree.
VCF or BAM format is supported for NRY haplogroup classification.
VCF files can be as input by the –-vcf option and only one VCF is required for analysis.
Multiple BAM files can be as input by adding -–bam option or providing a list of bam files.
'''


def classify_parser():

    parser = argparse.ArgumentParser('classify', description='(c) Y-LineageTracker: NRY haplogroup classification')
    # function used for classifying NRY haplogroups
    parser.add_argument('classify',
                        help='Classify NRY haplogroups and construct lineages.')
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
                        help='sample: A list containing all individuals used for NRY haplogroup classification')
    # optional, version of reference genome: 37 or 38
    parser.add_argument('-b', '--build',
                        required=True,
                        type=str,
                        action='store',
                        choices=['37', '38', 'Ref'],
                        help='build: The build version of the reference genome used for NRY haplogroup classification.')
    # optional, filter females by misisng rate
    parser.add_argument('--filter',
                        required=False,
                        type=str,
                        action='store',
                        nargs='?',
                        const='0.75',
                        help='filter: Filter female individuals by setting a cutoff of missing rate, the sample with a missing rate value greater than this value will be filtered. The default value is 0.75. If the optimal missing rate value is unknown, you can set the auto parameter to filter female individuals automatically.')
    # optional, check the data is whole genoeme sequence or SNP-array
    parser.add_argument('--chip',
                        required=False,
                        action='store_true',
                        help='chip: Declare the input data is SNP-array. This option is only used for the VCF file.')
    # optional, check whether classify haplogroup only from SNP data
    parser.add_argument('--snp-only',
                        required=False,
                        dest='snp',
                        action='store_true',
                        help='snp-only: Classify NRY haplogroups only from SNP data. This option is only used for the VCF file.')
    # optional, check whether output all matched mutation in classification
    parser.add_argument('--mut-info',
                        required=False,
                        dest='mutation',
                        action='store_true',
                        help='mutation: Output all matched haplogroup mutations of each individual.')
    # optional, check whether classify from approximate branch
    parser.add_argument('-a', '--approximate',
                        required=False,
                        action='store_true',
                        help='approximate: Classify NRY haplogroups which are in uncertain tree branches.')
    # optional, construct a tree from haplogroup results
    parser.add_argument('--phylo',
                        required=False,
                        type=str,
                        action='store',
                        nargs='?',
                        const='newick',
                        choices=['newick', 'nexus', 'phyloxml', 'nexml'],
                        help='phylo: Construct a phylogenetic tree based on NRY haplogroup classification results. This option is only used for the VCF file.')
    # optional, customization of mutation markers
    parser.add_argument('--ref',
                        required=False,
                        type=str,
                        action='store',
                        help='ref: The mutation markers provided by users.')
    # optional, the prefix of output
    parser.add_argument('-o', '--output',
                        required=False,
                        type=str,
                        action='store',
                        help='output: The prefix of output files.')

    args = parser.parse_args()

    return args


# print program information and write to log file
def set_log(input, output, format, cutoff, args_log):

    logger = logging.getLogger()
    logger.setLevel(level=logging.INFO)

    handler = logging.FileHandler(output, mode='w')
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter('[%(asctime)s] - [%(levelname)s]: %(message)s')
    handler.setFormatter(formatter)

    console = logging.StreamHandler()
    console.setLevel(logging.INFO)

    logger.addHandler(handler)
    logger.addHandler(console)

    log_info = ['[Y-LineageTracker] [Y-Haplogroup Classification]',
                '[Y-LineageTracker] Run Date: ' + time.asctime(time.localtime(time.time())),
                '[Y-LineageTracker] Input File: %s' % input,
                '[Y-LineageTracker] File Format: %s' % format,
                '[Y-LineageTracker] Reference Build: %s' % str(args_log.build),
                '[Y-LineageTracker] Data Type: WGS']

    if args_log.chip:
        log_info[4] = '[Y-LineageTracker] Data Type: SNP-array'
    if args_log.sample:
        log_info.append('[Y-LineageTracker] List of samples provided')
    else:
        if not cutoff is None:
            if cutoff:
                log_info.append('[Y-LineageTracker] Missing rate cutoff: %.2f' % cutoff)
            else:
                log_info.append('[Y-LineageTracker] Filter female individuals automatically')
    if args_log.mutation:
        log_info.append('[Y-LineageTracker] Output all matched mutation in classification')
    if args_log.phylo:
        log_info.append('[Y-LineageTracker] Output phylogenetic tree')

    print('\n')
    for i in log_info:
        logger.info(i)


# names of output files
def check_classify_output(mutation, output):

    from FilesIO import get_out_path
    path = get_out_path(output, 'Classification')
    haplogroup_file = path + '.hapresult.hg'
    lineage_file = path + '.lineageresult.txt'
    log_file = path + '.ClassificationLog.log'
    if mutation:
        mutation_file = path + '.mutationresult.txt'
        return haplogroup_file, lineage_file, mutation_file, log_file
    else:
        return haplogroup_file, lineage_file, log_file


def main():

    start = time.perf_counter()
    arguments = classify_parser()

    from FilesIO import CommonData, check_cutoff, check_ref_marker, check_BAM_VCF_input, check_samples, time_count
    # check input files
    input, format = check_BAM_VCF_input(arguments.vcf, arguments.bam, arguments.cram)

    # check cutoff: cutoff is a float between 0 and 1, else False
    cutoff = check_cutoff(arguments.filter, format)

    # set of output files
    output_set = check_classify_output(arguments.mutation, arguments.output)

    # check reference marker
    arguments.ref, arguments.build = check_ref_marker(arguments.ref, arguments.build)

    # set log file
    set_log(input, output_set[-1], format, cutoff, arguments)

    # check whether provide sample list
    sample = check_samples(arguments.sample)

    # read ISOGG haplogroup reference data
    common_data = CommonData()
    ref_info, ref_ind_info, approximate_haps = common_data.read_ref_info(arguments.snp, arguments.approximate, arguments.ref)

    # filter samples
    from ProcessData import ReadFilterData
    read_filter = ReadFilterData.ReadFilterData(input, sample)
    if format == 'bam':
        data_info, inds = read_filter.read_bam(ref_info, arguments.build, 'snp', format)
    elif format == 'cram':
        data_info, inds = read_filter.read_bam(ref_info, arguments.build, 'snp', format)
    else:
        data_info, inds = read_filter.read_filter_variant(cutoff, format)

    # match input data to ISOGG reference
    match = MatchVariant.MatchVariant(arguments.build)
    matched_info, header_cut = match.match_variant_snp(ref_info, data_info, format)

    # classify haplogroups and build lineage track
    tracker = TrackLineage.TrackLineage(matched_info, ref_ind_info, approximate_haps, arguments.build, format, arguments.chip)
    hap_data = tracker.classify_haplogroup(inds, matched_info, ref_info, output_set, arguments.mutation, header_cut)

    # construct phylogenetic tree
    if arguments.phylo:
        if arguments.vcf:
            if len(inds) > 1:
                from PhyloHaplogroup import check_phylo_output, PhyloHaplogroup
                output_set = check_phylo_output(arguments.phylo)
                phylo = PhyloHaplogroup(hap_data, output_set)
                hap_df, base, trunk_length, max_length, samples = phylo.build_base_tree()
                base = phylo.phylo_internal(input, format, base, samples)
                phylo.output_tree(base, arguments.phylo)
            else:
                logger = logging.getLogger()
                logger.warning('[Y-LineageTracker] Cannot construct phylognetic tree since there is only one individual')
        else:
            logger = logging.getLogger()
            logger.warning('[Y-LineageTracker] Constructing phylogenetic tree based on NRY haplogroup results only supports for VCF in current version')

    time_count(start)


if __name__ == '__main__':
    main()
