import os
import sys
import time
import logging
import argparse
import warnings
import itertools

import pandas as pd


class CommonData(object):
    '''
    This class is used to read common data for specific analysis
    '''

    def __init__(self, snp_only=False):

        warnings.filterwarnings('ignore')
        self.data_path = os.path.split(os.path.realpath(__file__))[0] + '/Data/'

    def read_region_info(self):

        region_info = pd.read_csv(self.data_path + 'Y_chromosome_region.bed',
                                  header=0,
                                  sep='\t',
                                  index_col=False,
                                  encoding='utf-8')

        return region_info

    def read_ref_info(self, snp_only, approximate, ref):

        if ref:
            import numpy as np
            from GetConfig import getConfig
            config = getConfig()
            common_trunk = config.get('HaplogroupTree', 'CommonTrunk').split(',')

            def get_approx(x):
                if x['Haplogroup'].endswith('~'):
                    if x['Haplogroup'].strip('~') in hgs:
                        return np.nan
                    else:
                        return 'Approximate'
                else:
                    return np.nan

            def get_approx_hap(x):
                if not x['Haplogroup'].endswith('~'):
                    return x['Haplogroup']
                else:
                    if pd.isna(x['ApproximateInfo']):
                        return x['Haplogroup']
                    else:
                        return x['Haplogroup'].strip('~')

            def get_main(x):
                if len(x['Haplogroup']) == 1:
                    return 'Main'
                elif x['Haplogroup'] in common_trunk:
                    return 'Main'
                else:
                    return np.nan

            # columns in custom data
            required_columns = ['Mutation', 'Haplogroup', 'Pos', 'Variant']
            # read custom data
            ref_info = pd.read_csv(ref,
                                   header=0,
                                   sep='\t',
                                   encoding='utf-8',
                                   dtype='object')
            if ref_info.columns.tolist() != required_columns:
                raise KeyError('Header of ref file not correct, please check the header is: Mutation, Haplogroup, Pos, Variant')

            # fill data
            ref_info['rs'] = '.'
            ref_info['BuildNa'] = '.'
            ref_info['MutationInfo'] = '.'
            ref_info['KeyInfo'] = '.'
            ref_info['ReferenceInfo'] = '.'
            hgs = sorted(set(ref_info[ref_info['Haplogroup'].map(lambda x: not x.endswith('~'))]['Haplogroup'].tolist()))
            ref_info['ApproximateInfo'] = ref_info.apply(get_approx, axis=1)
            ref_info['Haplogroup'] = ref_info.apply(get_approx_hap, axis=1)
            ref_info['MainInfo'] = ref_info.apply(get_main, axis=1)

            ref_info.rename(columns={'Pos': 'BuildRef'}, inplace=True)

            ref_info = ref_info[['Mutation', 'Haplogroup', 'rs', 'BuildRef', 'BuildNa', 'Variant',
                                 'MutationInfo', 'KeyInfo', 'ReferenceInfo', 'ApproximateInfo', 'MainInfo']]
        else:
            ref_info = pd.read_csv(self.data_path + 'ISOGG_haplogroup_info.csv',
                                   header=0,
                                   sep=',',
                                   encoding='utf-8',
                                   dtype='object')
        if not approximate:
            ref_info = ref_info[ref_info['Haplogroup'].map(lambda x: not x.endswith('~'))]
        approximate_haps = sorted(set(ref_info[ref_info['ApproximateInfo']=='Approximate']['Haplogroup']))
        if snp_only:
            ref_info = ref_info[ref_info['Variant'].map(lambda x: len(x) == 1)]
        ref_ind_info = ref_info[ref_info['ReferenceInfo']=='Ref']

        return ref_info, ref_ind_info, approximate_haps


    def read_pop_info(self):

        pop_info = pd.read_csv(self.data_path + 'ISOGG_population_index.csv',
                               header=0,
                               sep=',',
                               encoding='utf-8')

        return pop_info

    def read_str_info(self):

        str_info = pd.read_csv(self.data_path + 'Y_STR_panel.csv',
                               header=0,
                               sep=',',
                               encoding='utf-8')

        return str_info

    def read_str_prediction_panel(self):

        Y_STR_precidtion_panel = pd.read_csv(self.data_path + 'Y_STR_precidtion_panel.csv',
                                             header=0,
                                             sep=',',
                                             encoding='utf-8')

        return Y_STR_precidtion_panel

    def read_prior_calibration(self):

        prior_calibration = pd.read_csv(self.data_path + 'haplogroup_calibration.csv',
                                        header=0,
                                        index_col=0,
                                        sep=',',
                                        encoding='utf-8')

        return prior_calibration

    def read_tree_info(self):

        tree_file = open(self.data_path + 'haplogroup_tree_2018.txt', 'r')
        tree_info = list(map(lambda x: x.replace('\n', ''), tree_file.readlines()))
        tree_file.close()

        return tree_info


# get output path
def get_out_path(output, name):

    if output:
        if not os.path.isdir(output):
            if os.path.isdir(os.path.split(output)[0]) or os.path.split(output)[0] == '':
                path = output
            else:
                raise FileNotFoundError('No such directory: %s' % os.path.split(output)[0])
        else:
            raise FileNotFoundError('%s is directory' % output)
    else:
        path = os.getcwd() + '/' + name

    return path


# get STR panel
def get_str_from_panel(panel):

    common_data = CommonData()
    str_info = common_data.read_str_info()
    # choose which kit will be used to genotype Y-STR
    if panel != 'all':
        str_info = str_info[str_info['STR'].map(lambda x: x != '.')]
    if panel != 'named':
        str_info = str_info[str_info['Kit'].map(lambda x: not pd.isna(x))]
        if panel == 'minimal':
            str_info = str_info[str_info['Kit'].map(lambda x: 'Minimal' in x.split(', '))]
        elif panel == 'ppy':
            str_info = str_info[str_info['Kit'].map(lambda x: 'PowerPlex_Y' in x.split(', '))]
        elif panel == 'pp23':
            str_info = str_info[str_info['Kit'].map(lambda x: 'PowerPlex_Y23' in x.split(', '))]
        elif panel == 'yf':
            str_info = str_info[str_info['Kit'].map(lambda x: 'Yfiler' in x.split(', '))]

    return str_info


# check whether set cutoff of filter
def check_cutoff(filter, format):

    if format == 'bam':
        if filter:
            logger = logging.getLogger()
            logger.warning('[Y-LineageTracker] [Warning] filter option only used for VCF' % len(overlap))
        cutoff = None
    else:
        if filter == None:
            cutoff = None
        elif filter == 'auto':
            cutoff = False
        else:
            try:
                cutoff = float(filter)
            except:
                raise argparse.ArgumentTypeError('Cutoff of missing rate should be a float between 0 and 1')
            if cutoff < 0 or cutoff > 1:
                raise argparse.ArgumentTypeError('Cutoff of missing rate should be a float between 0 and 1')

    return cutoff


# check reference marker file
def check_ref_marker(ref, build):

    if not ref:
        if build == 'Ref':
            raise argparse.ArgumentTypeError('No reference marker file specified')
    else:
        if build != 'Ref':
            logger = logging.getLogger()
            logger.warning('[Y-LineageTracker] [Warning] Not specify build argument Ref, use %d' % build)

    return ref, build

# check bam or vcf input
def check_BAM_VCF_input(vcf, bam, cram):

    if vcf:
        input = vcf
        compressed = check_compress(input)
        if compressed:
            format = 'gzvcf'
        else:
            format = 'vcf'
    else:
        if bam:
            input = ', '.join(bam)
            format = 'bam'
        else:
            input = ', '.join(cram)
            format = 'cram'


    return input, format


# check sequence alignment file or matrix file
def check_seq_matrix_input(seq, pre_format, matrix):

    if seq:
        input = seq
        if pre_format:
            format = pre_format
        else:
            raise argparse.ArgumentTypeError('Fomrat should be specified for sequence file by --seq-format')
    else:
        input = matrix
        format = 'matrix'

    return input, format


# check whther file is compressed
def check_compress(input_file):

    if os.path.isfile(input_file):
        try:
            open(input_file).read()
            compressed = False
        except UnicodeDecodeError:
            compressed = True
    else:
        raise FileNotFoundError('%s not found, please check your file name' % input_file)
        sys.exit()

    return compressed


# check whether sample file is provided
def check_samples(male_file):

    if male_file:
        try:
            with open(male_file, 'r') as f:
                males = f.read().splitlines()
        except:
            raise FileNotFoundError('Unable to open male file, please check your file name')
    else:
        males = []

    return males


# check header and value type of input data
def check_hap_input(input_file, type, format='matrix'):

    if type == 'haplotype' and format != 'matrix':
            from ProcessData import ConvertData
            convert = ConvertData.ConvertData(input_file)
            input_data = convert.data_convert(format, 'matrix')
    else:
        input_data = pd.read_csv(input_file,
                                 sep='\t',
                                 header=0,
                                 encoding='utf-8')
    input_data = input_data.apply(lambda x: x.str.strip() if x.dtype == "object" else x) # strip whitespace

    if type == 'haplogroup':
        if input_data.columns.size != 2 or input_data.columns.tolist() != ['SampleID', 'Haplogroup']:
            raise KeyError('Header of haplogroup file not correct, please check the header is: SampleID, Haplogroup')
        else:
            input_data['Haplogroup'] = input_data['Haplogroup'].map(lambda x: x.strip('~'))
    elif type == 'haplotype':
        if input_data.columns.size < 2 or input_data.columns.tolist()[0] != 'SampleID':
            raise KeyError('Header of haplotype file not correct, please check the header is: SampleID, Haplotype Names')

    return input_data


# check population file
def check_population(population_file):

    try:
        population_data = pd.read_csv(population_file,
                                      sep='\t',
                                      header=0,
                                      encoding='utf-8')
        population_data = population_data.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    except:
        raise FileNotFoundError('%s not found, please check your file name' % population_file)
    condition_a = population_data.columns.size == 2 and population_data.columns[0:2].tolist() == ['SampleID', 'Population']
    condition_b = population_data.columns.size == 3 and population_data.columns[0:3].tolist() == ['SampleID', 'Population', 'Group']
    condition_c = population_data.columns.size == 3 and population_data.columns[0:3].tolist() == ['SampleID', 'Population', 'Color']
    condition_d = population_data.columns.size == 3 and population_data.columns[0:3].tolist() == ['SampleID', 'Population', 'Shape']
    condition_e = population_data.columns.size == 4 and population_data.columns[0:4].tolist() == ['SampleID', 'Population', 'Group', 'Color']
    condition_f = population_data.columns.size == 4 and population_data.columns[0:4].tolist() == ['SampleID', 'Population', 'Group', 'Shape']
    condition_g = population_data.columns.size == 4 and population_data.columns[0:4].tolist() == ['SampleID', 'Population', 'Color', 'Shape']
    condition_h = population_data.columns.size == 5 and population_data.columns[0:5].tolist() == ['SampleID', 'Population', 'Group', 'Color', 'Shape']
    if not (condition_a or condition_b or condition_c or condition_d or condition_e or condition_f or condition_g or condition_h):
        raise KeyError('Header of population file not correct, please check the header is: SampleID, Population or SampleID, Population, Group')

    return population_data


# check population header
def check_pop_header(population_data):

    if 'Group' in population_data.columns.tolist():
        population_data = population_data[['SampleID', 'Population', 'Group']]
    else:
        population_data = population_data[['SampleID', 'Population']]

    return population_data


# set SampleID as index of haplotype data
def set_haplotype_index(haplotype_data, format):

    from GetConfig import getConfig
    config = getConfig()
    skip_letters = config.get('Statistics', 'SkipLetters').split(',')

    haplotype_data = haplotype_data.set_index(['SampleID'], drop=True)
    if format == 'matrix':
        if haplotype_data.applymap(lambda x: isinstance(x, int) or x.isdigit() or x in skip_letters).sum().sum() != haplotype_data.size:
            raise KeyError('Unexpected value in input haplotype matrix data')

    return haplotype_data


# check whther individuals are overlapped in two files
def check_overlap(sample_data, population_data):

    logger = logging.getLogger()

    # sort and drop duplicated individuals
    sample_data = sample_data.sort_values(by='SampleID').drop_duplicates().reset_index(drop=True)
    population_data = population_data.sort_values(by='SampleID').drop_duplicates().reset_index(drop=True)
    # whther samples are same of two files
    if sample_data.index.size == population_data.index.size:
        if all(sample_data['SampleID']==population_data['SampleID']):
            pass
        else:
            # only keep matched indvisuals
            TF = sample_data['SampleID']==population_data['SampleID']
            sample_data = sample_data[TF]
            population_data = population_data[TF]
            overlap_num = sum(sample_data['SampleID']==population_data['SampleID'])
            logger.info('[Y-LineageTracker] [Warning] Only %d individuals overlaped in haplogroup file and population file' % len(overlap_num))
    else:
        overlap = sorted(list(set(sample_data['SampleID'])&set(population_data['SampleID'])))
        if len(overlap) > 0:
            # only keep matched indvisuals
            sample_data = sample_data[sample_data['SampleID'].map(lambda x: x in overlap)]
            population_data = population_data[population_data['SampleID'].map(lambda x: x in overlap)]
            logger.info('[Y-LineageTracker] [Warning] Only %d individuals overlaped in haplogroup file and population file' % len(overlap))
        else:
            # No matched indvisuals
            logger.error('[Y-LineageTracker] Format error or There are no overlaped individuals in haplogroup/haplotype file and population file')
            sys.exit()
    population_data = population_data.sort_values(by=population_data.columns.tolist()[::-1][1:]).reset_index(drop=True)

    return sample_data, population_data


# set default colors
def set_color_num(num, get_cmap=False):

    from matplotlib import cm
    from GetConfig import getConfig

    config = getConfig()
    color_num = num
    mim_color_num = config.getint('PlotParameters', 'MinColorNum')
    median_color_num = config.getint('PlotParameters', 'MedianColorNum')
    max_color_num = config.getint('PlotParameters', 'MaxColorNum')
    if color_num <= mim_color_num:
        cmap = cm.get_cmap('Set1').colors
    elif color_num > mim_color_num and color_num <= median_color_num:
        cmap = cm.get_cmap('tab20').colors
    elif color_num > median_color_num and color_num <= max_color_num:
        cmap = cm.get_cmap('tab20').colors + cm.get_cmap('tab20b').colors
    colors = list(cmap)[:color_num]

    if get_cmap:
        return colors, cmap
    else:
        return colors


# haplogroup normalization
def haplogroup_normalize(haplogroup):

    upper_trunk = config.get('HaplogroupTree', 'UpperTrunk')

    if not haplogroup[0].isupper():
        haplogroup = haplogroup[0].upper()+haplogroup[1:]
    if not len(haplogroup) == 1:
        if haplogroup.upper() in upper_trunk:
            haplogroup = haplogroup.upper()
        else:
            if not haplogroup[1:].islower():
                haplogroup = haplogroup[0] + haplogroup[1:].lower()

    return haplogroup


# get end time to calculate total run time
def time_count(start, mid=False):

    logger = logging.getLogger()

    interval = time.perf_counter() - start
    if interval > 60:
        minute = int(interval/60)
        second = interval - minute * 60

    if mid:
        if interval > 60:
            logger.info('[Y-LineageTracker] Time: %dmin%.2fs' % (minute, second))
        else:
            logger.info('[Y-LineageTracker] Time: %.2fs' % (interval))
    else:
        if interval > 60:
            logger.info('[Y-LineageTracker] Total run time: %dmin%.2fs' % (minute, second))
        else:
            logger.info('[Y-LineageTracker] Total run time: %.2fs' % (interval))
