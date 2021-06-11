import os
import sys
import time
import logging
import argparse
import itertools
import matplotlib; matplotlib.use('Agg')
import multiprocessing
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from operator import add
from functools import partial
from collections import Counter
from matplotlib import font_manager
from GetConfig import getConfig


config = getConfig()


'''
Description:
This module is used for network analysis.
Network analysis is an efficient method to show relationship of haplotypes.
Y haplotype data is required for analysis.
It is recommended to use Y-STR haplotype under a specific haplogroup.
Network analysis can show the evolutionary relationship between haplotypes.
This command can output network figure, and fdi file.
'''


def net_parser():

    parser = argparse.ArgumentParser('net', description='(c) Y-LineageTracker: Network analysis')
    # function used for network analysis
    parser.add_argument('net',
                        help='Perform network analysis from Y haplotype data.')
    # required, input file, alignment sequence or matrix
    input = parser.add_mutually_exclusive_group(required=True)
    input.add_argument('--seq',
                        type=str,
                        action='store',
                        help='seq: Haplotype data in sequence alignment format')
    input.add_argument('--matrix',
                        type=str,
                        action='store',
                        help='martix: Y-STR haplotype data in matrix format')
    # optional, format of sequence alignment
    parser.add_argument('--seq-format',
                        required=False,
                        type=str,
                        dest='format',
                        action='store',
                        choices=['fasta', 'phylip', 'nexus', 'meg', 'vcf'],
                        help='--seq-format: The format of sequence file. This option is only required for the sequence alignment file. ')
    # optional, population file
    parser.add_argument('-p', '--population',
                        required=False,
                        type=str,
                        action='store',
                        help='population: A file containing sample ID and population information of each individual.')
    # optional, type of tree
    parser.add_argument('--tree-type',
                        required=False,
                        type=str,
                        dest='treetype',
                        action='store',
                        default='mjn',
                        choices=['mjn', 'msn'],
                        help='tree-type: The algorithm used to construct a network tree.')
    # optional, how to code missing site
    parser.add_argument('--gap',
                        required=False,
                        type=str,
                        action='store',
                        default='allele',
                        choices=['allele', 'missing'],
                        help='gap: The gap coding method. You can choose to code gaps as a new type of allele or missing data')
    # optional, the search for new median nodes
    parser.add_argument('--search',
                        required=False,
                        type=int,
                        action='store',
                        default=0,
                        help='search: How the search for new median nodes is performed: the larger this parameter, the wider the search. ')
    # optional, output fdi file
    parser.add_argument('--fdi',
                        required=False,
                        action='store_true',
                        help='fdi: Output fdi file that can be used as the input file in Network software, to show network plot or estimate TMRCA directly.')
    # optional, filter sites by misisng rate value
    parser.add_argument('--filter',
                        required=False,
                        type=float,
                        action='store',
                        default=0.7,
                        help='filter: The cutoff value of the site to be filtered, sites with a missing rate value greater than this value will be filtered')
    # optional, how to show alternative link in figure
    parser.add_argument('--alt',
                        required=False,
                        type=str,
                        action='store',
                        default='hide',
                        choices=['hide', 'show', 'dashed'],
                        help='alt: Whether to show alternative links between nodes in the network figure.')
    # optional, size of median nodes
    parser.add_argument('--mediansize',
                        required=False,
                        type=float,
                        action='store',
                        default=0.0,
                        help='mediansize: The median node size in the median-joining network compared to base size in the network figure. The default is 0.')
    # optional, keep all nodes in same size
    parser.add_argument('--euqal-size',
                        required=False,
                        dest='euqalsize',
                        action='store_true',
                        help='euqal-size: Keep all nodes the same size in the network figure.')
    # optional, the prefix of output
    parser.add_argument('-o', '--output',
                        required=False,
                        type=str,
                        action='store',
                        help='output: The prefix of output files. ')

    args = parser.parse_args()

    return args


# check missing rate is filter argument is used
def check_misisng_rate(missing_rate):

    if missing_rate < 0 or missing_rate> 1:
        print('[Y-LineageTracker] [Error] Cutoff of missing rate should be a float between 0 and 1')
        sys.exit()


# print program information and write to log file
def set_log(args_log, log_file):

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

    log_info = ['[Y-LineageTracker] [Network]',
                '[Y-LineageTracker] Run Date: ' + time.asctime(time.localtime(time.time())),
                '[Y-LineageTracker] Haplotype File: %s' % args_log.seq,
                '[Y-LineageTracker] Haplotype File Format: %s' % args_log.format,
                '[Y-LineageTracker] Gap euqal to %s' % args_log.gap,
                '[Y-LineageTracker] Searching parameter: %d' % args_log.search,
                '[Y-LineageTracker] Missing rate of allele to be filtered: %s' % args_log.filter,
                '[Y-LineageTracker] Alternative link type : %s' % args_log.alt,
                '[Y-LineageTracker] Median node size: %s' % args_log.mediansize]

    if args_log.euqalsize:
        log_info.append('[Y-LineageTracker] Keep all node the same size')
    else:
        log_info.append('[Y-LineageTracker] Plot node according to sample size')

    print('\n')
    for i in log_info:
        logger.info(i)


class HaploNetwork(object):
    '''
    This class is used for network analysis.
    Main functions of this class include:
    1. construct median-joining network from haplotypes
    2. plot network by networkx package
    3. write network to fdi file
    4. output network info
    '''

    def __init__(self, gap, path):

        self.logger = logging.getLogger()
        self.skip_letters = config.get('Statistics', 'SkipLetters').split(',')
        self.gap = gap
        self.path = path
        self.net_info = open(self.path+'.info.net', 'w')

    # convert STR numbers to 0-1 code haplotype
    def _convert_allel_to_binary(self, haplotype_data, filter):

        # create a empty dataframe to save haplotypes
        sequence_df = pd.DataFrame(index=haplotype_data.index,
                                   columns=['Haplotype'])
        sequence_df['Haplotype'] = ''

        STR_code_pos = {}
        start_pos = 0
        hap_length = 0
        allel_num = haplotype_data.index.size
        for STR in haplotype_data.columns.tolist():
            str_allels = haplotype_data[STR].astype(str).tolist()
            # filter STR according to filter parameter
            if sum([str_allels.count(i) for i in self.skip_letters])/allel_num >= filter:
                continue
            if len(set(set(str_allels)&set(self.skip_letters))) > 1:
                str_allels = ['.' if i in self.skip_letters else i for i in str_allels]
            repeats_type = sorted(list(set(str_allels)))
            # convert STR repeat numbers to 0-1 code
            repeats_type_num = len(repeats_type)
            max_length = len('{0:b}'.format(repeats_type_num-1))
            binary_set = [] # creat a empty list to save 0-1 code
            for i in range(repeats_type_num):
                binary_num = '{0:b}'.format(i)
                binary_num_length = len(binary_num)
                if binary_num_length < max_length:
                    binary_num = (max_length-binary_num_length)*'0'+binary_num
                binary_set.append(binary_num)
            if self.gap == 'missing':
                replaced_allels = ['-'*max_length if i == '.' else binary_set[repeats_type.index(i)] for i in str_allels]
            else:
                replaced_allels = [binary_set[repeats_type.index(i)] for i in str_allels]
            sequence_df['Haplotype'] = list(map(add, sequence_df['Haplotype'].tolist(), replaced_allels))
            for pos in range(start_pos, start_pos+max_length):
                STR_code_pos[pos] = STR
            start_pos += max_length
            hap_length += 1

        self.net_info.write('#HapLength:%d\n\n' % hap_length)

        return sequence_df, STR_code_pos

    # calculate frequency of all kinds of haplotypes and the proportion in each population
    def _calculate_distance_matrix(self, sequence_df):

        # get haplotypes and name for them
        haplotypes = sorted(list(set(sequence_df['Haplotype'])))
        sp_sequence_df = pd.DataFrame(index=range(len(haplotypes)), columns=['Haplotype'])
        sp_sequence_df['Haplotype'] = haplotypes

        # calculate distance between two haplotypes
        distance_matrix = pd.DataFrame(index=range(len(haplotypes)), columns=range(len(haplotypes)))
        hap_num = 0
        for hap1_name, hap1 in enumerate(haplotypes):
            hap_num += 1
            for hap2_name, hap2 in enumerate(haplotypes[hap_num:]):
                if self.gap == 'missing':
                    distance = sum(l1 != l2 for l1, l2 in zip(hap1, hap2) if l1 not in self.skip_letters and l2 not in self.skip_letters)
                    if distance == 0:
                        distance = 1
                else:
                    distance = sum(l1 != l2 for l1, l2 in zip(hap1, hap2))
                hap2_name = hap2_name + hap_num
                distance_matrix.at[hap2_name, hap1_name] = distance

        self.net_info.write('#HapSamples:\n')

        sample_sp_code = {}
        for i in haplotypes:
            idx = sequence_df['Haplotype'].tolist().index(i)
            sample = sequence_df.index[idx]
            sample_sp_code[sample] = i
            all_samples = ','.join(sequence_df[sequence_df['Haplotype']==i].index.tolist())
            self.net_info.write('%d\t%s\n' %(idx, all_samples))

        self.net_info.write('\n')

        return distance_matrix, sp_sequence_df, sample_sp_code

    # get msn tree
    def _calculate_minimum_spanding_tree(self, n, distance_matrix, return_pruned=True):

        if n < 2:
            self.logger.error('[Y-LineageTracker] Program stopped since there are less than 2 types of haplotypes')
            sys.exit()

        # create a empty edge matrix to save hapotype and their distance
        edge_num = n-1
        edge_matrix = np.empty(shape=(edge_num, 3))
        edge_matrix[:] = np.nan

        # get the distance from small to large between of two haplotypes
        forest = np.linspace(0, n-1, n)

        add_num = 0
        while add_num < edge_num:
            ci, ri = np.unravel_index(np.nanargmin(distance_matrix.T.values), distance_matrix.T.shape)
            row_record = forest[ri]
            col_record = forest[ci]
            if row_record != col_record:
                edge_matrix[add_num] = [ci, ri, distance_matrix.at[ri, ci]] # col hap, row hap, step
                forest[forest==col_record] = row_record
                add_num += 1
            distance_matrix.at[ri, ci] = np.nan

        if return_pruned:
            return edge_matrix, distance_matrix
        else:
            return edge_matrix

    # get median haplotypes
    @staticmethod
    def calculate_median_vector(triplet, sp_sequence_df, max_n_cost, epsilon, n, gap, new_link_data, shared_median_haps, initial_lambda, n_cost, lock):

        median_haps = pd.DataFrame(columns=['Haplotype'])
        temp_link_data = np.empty(shape=(0, 3))
        median_num = n

        for item in triplet:
            hap_in_order = item
            # check whether two haplotypes can be linked
            condition1 = any(hap_in_order[0]==np.array(new_link_data)[:, 0]) and any(hap_in_order[1]==np.array(new_link_data)[:, 1])
            condition2 = any(hap_in_order[0]==np.array(new_link_data)[:, 0]) and any(hap_in_order[2]==np.array(new_link_data)[:, 1])
            condition3 = any(hap_in_order[1]==np.array(new_link_data)[:, 0]) and any(hap_in_order[2]==np.array(new_link_data)[:, 1])
            if (condition1+condition2+condition3) < 2:
                continue
            # get alleles that are variable
            sub_seq = sp_sequence_df.loc[list(hap_in_order)]
            # get the median haplotype
            sub_seq_length = len(sub_seq.at[hap_in_order[0], 'Haplotype'])
            median_seq = ''.join([Counter(sub_seq['Haplotype'].map(lambda x: x[i])).most_common(1)[0][0] for i in range(sub_seq_length)])
            # check whether median sequence exist
            if median_seq in sp_sequence_df['Haplotype'].tolist():
                continue
            if not median_haps.empty:
                if median_seq in shared_median_haps:
                    continue
            # calculate the distance between median sequence and haplotypes
            if gap == 'missing':
                median_distance = [sum (l1 != l2 for l1, l2 in zip(i, median_seq) if l1 not in self.skip_letters and l2 not in self.skip_letters) for i in sub_seq['Haplotype']]
                if 0 in median_distance:
                    median_distance = [1 if dist==0 else dist for dist in median_distance]
            else:
                median_distance = [sum (l1 != l2 for l1, l2 in zip(i, median_seq)) for i in sub_seq['Haplotype']]
            connection_cost = sum(median_distance)
            # update n cost
            lock.acquire()
            n_cost.value += 1
            if n_cost.value > max_n_cost:
                print('[Y-LineageTracker] Program stopped since there are too many connnections')
                sys.exit()
            if connection_cost < initial_lambda.value:
                # update initial lambda
                initial_lambda.value = connection_cost
            if connection_cost <= initial_lambda.value + epsilon:
                hap_to_add1 = hap_in_order
                hap_to_add2 = [median_num]*3
                step_to_add = median_distance
                link_data_to_add = np.array([hap_to_add1, hap_to_add2, step_to_add]).T
                temp_link_data = np.concatenate([temp_link_data, link_data_to_add])
                median_haps.loc[median_num] = median_seq
                median_num += 1
                if median_seq not in shared_median_haps:
                    shared_median_haps.append(median_seq)
            lock.release()

        return median_haps, temp_link_data

    # construct network from distance matrix
    def _calculate_network(self, distance_matrix, sp_sequence_df, tree_type, epsilon):

        # initialization
        n = sp_sequence_df.index.size
        haplotype_length = len(sp_sequence_df.at[0, 'Haplotype'])
        median_num = n
        MST, pruned_distance_matrix = self._calculate_minimum_spanding_tree(n ,distance_matrix)

        if n < 3 or tree_type == 'msn':
            self.logger.warning('There are less than 3 types of haplotypes')
            alt_data = np.empty(shape=(0, 3))
            return MST, alt_data

        MST_step = MST[:, 2]
        lower_limit = min(MST_step)
        upper_limit = max(MST_step) + epsilon
        new_link_data = []
        for i in pruned_distance_matrix.columns:
            for j in pruned_distance_matrix.index:
                step = pruned_distance_matrix.at[j, i]
                if (not np.isnan(step)) and step > lower_limit and step <= upper_limit:
                    new_link_data.append([i, j, step])
        if new_link_data:
            new_link_data = np.concatenate([MST, np.array(new_link_data)])
        else:
            new_link_data = MST

        start_n_cost = 0
        start_initial_lambda = float('inf')
        tripleted = []

        while True:

            # get triplet for calculation
            triplet_new = list(itertools.combinations(range(n), 3))
            triplet = [i for i in triplet_new if i not in tripleted]
            tripleted.extend(triplet)

            process_num = 10
            pool = multiprocessing.Pool(processes=process_num)

            # global variables to be shared
            manager = multiprocessing.Manager()
            n_cost = manager.Value('i', start_n_cost)
            initial_lambda  = manager.Value('f', start_initial_lambda)
            shared_new_link_data = manager.list(new_link_data)
            shared_median_haps = manager.list([])
            LOCK = manager.Lock()

            # split triplet to sections
            triplet_interval = np.linspace(0, len(triplet), process_num, dtype=int)
            #triplet_sections = [triplet[triplet_interval[n]: triplet_interval[n+1]] if n == 0 else triplet[(triplet_interval[n]+1): triplet_interval[n+1]] for n in range(process_num-1)]
            triplet_sections = [triplet[triplet_interval[n]: triplet_interval[n+1]] for n in range(process_num-1)]

            # multiprocess triplet data
            pt = partial(self.calculate_median_vector, # function
                         sp_sequence_df=sp_sequence_df, # variable for read
                         max_n_cost=40000, # variable for read
                         epsilon=epsilon, # variable for read
                         n=n, # variable for read
                         gap=self.gap, # variable for read
                         new_link_data=new_link_data, # variable for writing
                         shared_median_haps=shared_median_haps, # variable for writing
                         initial_lambda=initial_lambda, # variable for writing
                         n_cost=n_cost, # variable for writing
                         lock=LOCK) # process lock
            median_haps_list, temp_link_data_list = zip(*pool.map(pt, triplet_sections))

            pool.close()
            pool.join()

            # intergrate data
            median_haps = pd.concat(median_haps_list)
            temp_link_data = np.concatenate(temp_link_data_list)
            if median_haps.empty:
                self.logger.info('[Y-LineageTracker] Get all median haplotypes')
                break
            else:
                # code for median vectors
                start_num = sp_sequence_df.index.size
                end_num = start_num + median_haps.index.size
                median_haps.index = range(start_num, end_num)
                temp_link_data[:, 2] = [val for val in range(start_num, end_num) for _ in (0, 1, 2)]
                # add to haplotype information
                sp_sequence_df = pd.concat([sp_sequence_df, median_haps])
                new_link_data = np.concatenate([new_link_data, temp_link_data])
                # update for iteration
                start_n_cost = n_cost.value
                start_initial_lambda = initial_lambda.value

        sp_sequence_df = sp_sequence_df.drop_duplicates().reset_index()

        self.logger.info('[Y-LineageTracker] Calculating distance matrix...')
        # calculate and add median haplotype to distance matrix
        hap_num = 0
        for hap1_name, hap1 in zip(sp_sequence_df.index, sp_sequence_df['Haplotype']):
            hap_num += 1
            for hap2_name, hap2 in zip(sp_sequence_df.index[hap_num:], sp_sequence_df['Haplotype'][hap_num:]):
                if self.gap == 'missing':
                    distance = sum(l1 != l2 for l1, l2 in zip(hap1, hap2) if l1 not in self.skip_letters and l2 not in self.skip_letters)
                    if distance == 0:
                        distance = 1
                else:
                    distance = sum(l1 != l2 for l1, l2 in zip(hap1, hap2))
                distance_matrix.at[hap2_name, hap1_name] = distance


        N = sp_sequence_df.index.size
        dnet = distance_matrix.copy().applymap(lambda x: 0 if not np.isnan(x) else x)
        forest = np.linspace(0, N-1, N)
        network_data = np.empty(shape=(0, 3))
        alt_data = np.empty(shape=(0, 3))
        steps = sorted(list(set(filter(lambda x: x==x, distance_matrix.to_numpy().reshape(-1)))))

        edge_num = 0
        add_new_data = False

        self.logger.info('[Y-LineageTracker] Constructing relationship between nodes')
        while len(set(forest)) > 1:

            step_num = steps[edge_num]
            step_df = distance_matrix[distance_matrix==step_num]
            step_count = step_df.notnull().sum().sum()
            for i in range(step_count):
                ci, ri = np.unravel_index(np.nanargmin(step_df.T.values), step_df.T.shape)
                step_df.at[ri, ci] = np.nan
                row_record = forest[ri]
                col_record = forest[ci]
                if row_record != col_record:

                    for col in np.where(forest==col_record)[0]:
                        if ci > col:
                            num_col = dnet.at[ci, col]
                        elif ci < col:
                            num_col = dnet.at[col, ci]
                        else:
                            num_col = 0
                        for row in np.where(forest==row_record)[0]:
                            if ri > row:
                                num_row = dnet.at[ri, row]
                            elif ri < row:
                                num_row = dnet.at[row, ri]
                            else:
                                num_row = 0
                            dnet_num = num_col + num_row + step_num
                            dnet.at[row, col] = dnet_num
                    forest[forest==col_record] = row_record
                    add_new_data = True
                else:
                    if step_num < dnet.at[ri, ci]:
                        add_new_data = True
                        # update dnet
                        cols = np.append(network_data[network_data[:, 1]==ci, 0], network_data[network_data[:, 0]==ci, 1])
                        rows = np.append(network_data[network_data[:, 1]==ri, 0], network_data[network_data[:, 0]==ri, 1])
                        for col in np.append(ci, cols).astype(int):
                            if ci > col:
                                num_col = dnet.at[ci, col]
                            elif ci < col:
                                num_col = dnet.at[col, ci]
                            else:
                                num_col = 0
                            for row in np.append(ri, rows).astype(int):
                                if col == row:
                                    continue
                                else:
                                    if ri > row:
                                        num_row = dnet.at[ri, row]
                                    elif ri < row:
                                        num_row = dnet.at[row, ri]
                                    else:
                                        num_row = 0
                                    dnet_num = num_col + num_row + step_num
                                    if row > col:
                                        if dnet.at[row, col] > dnet_num:
                                            dnet.at[row, col] = dnet_num
                                    elif row < col:
                                        if dnet.at[col, row] > dnet_num:
                                            dnet.at[col, row] = dnet_num
                if add_new_data:
                    network_to_add = np.array([[ci, ri, step_num]])
                    network_data = np.concatenate([network_data, network_to_add])
                    add_new_data = False
            edge_num += 1
            if edge_num > 5000:
                self.logger.error('Program stopped since there are too many connnections')
                sys.exit()

        network_data = network_data.astype(int)
        net_MST = self._calculate_minimum_spanding_tree(N ,distance_matrix, False)
        net_MST = net_MST.astype(int)

        if len(network_data) > len(net_MST):
            alt_data = np.array([i for i in network_data if i[:2].tolist() not in net_MST[:, :2].tolist()])
            network_data = np.array([i for i in network_data if i[:2].tolist() in net_MST[:, :2].tolist()])
        else:
            alt_data = np.array([])

        self.logger.info('[Y-LineageTracker] Calculating median haplotypes finished')

        return network_data, alt_data, sp_sequence_df

    # prune network for visualization
    def _get_network_for_plot(self, network_data, alt_data, hap_num, alt_type):

        # specify network data according to type of alternative link
        if alt_data.size != 0:
            alt_type = alt_type
            pruned_network_data = network_data.copy()
            max_number = np.max(pruned_network_data[:, :2]) + 1
            while True:
                to_be_deleted_row = []
                for i in range(hap_num, max_number):
                    if np.count_nonzero(pruned_network_data[:, :2]==i) == 1:
                        row_index = np.where(pruned_network_data[:, :2]==i)[0]
                        to_be_deleted_row.append(row_index)
                pruned_network_data = np.delete(pruned_network_data, to_be_deleted_row, axis=0)
                if len(to_be_deleted_row) == 0:
                    break
            # get network data for visulization
            if alt_type == 'hide':
                network_data = pruned_network_data
        else:
            alt_type = 'hide'
            pruned_network_data = network_data.copy()

        self.net_info.write('#HapLinks:\n')
        [self.net_info.write('\t'.join(i)+'\n') for i in pruned_network_data.astype(str)]
        self.net_info.close()
        node_num = len(np.unique(network_data[:, :2]))

        return network_data, node_num

    # get mutation number
    def _get_mutations(self, network_line, sp_sequence_df, STR_code_pos):

        seq1 = sp_sequence_df.at[network_line[0], 'Haplotype']
        seq2 = sp_sequence_df.at[network_line[1], 'Haplotype']
        mutation = [STR_code_pos[n] for n, s in enumerate(zip(seq1, seq2)) if s[0] != s[1]]

        return mutation

    # write and output fdi file
    def _get_fdi_file(self, network_data, torso_num, node_num, sequence_df, sp_sequence_df, sample_sp_code, pos, hap_size, hap_freq, STR_code_pos, pops):

        fdi_file = open(self.path+'.fdi', 'w')
        medians = sorted(list(set(set(range(len(sample_sp_code), sp_sequence_df.index.size)) & set(pos.keys()))))

        # head info
        center_x,  center_y = np.mean(np.array(list(pos.values())).T, axis=1)
        fdi_file.write('SHOW_NODENAMES;TRUE\n'\
                       'NODES_PROPORTIONAL;TRUE\n'\
                       'SHOW_ONLY_TORSO;FALSE\n'\
                       'SHOW_MEDIANS;FALSE\n'\
                       'SHOW_MEDIAN_NAMES;TRUE\n'\
                       'SHOW_CHARACTERS;FALSE\n'\
                       'SHOW_LINES;FALSE\n'\
                       'MIN_CIRC_RADIUS;4\n'\
                       'MAX_CIRC_RADIUS;50\n'\
                       'MED_RADIUS;3\n'\
                       '\n'\
                       'SPECIAL_FILE_TYPE;FALSE\n'\
                       'DIAGRAM_CENTER;%d;%d\n'\
                       'IMAGE_CENTER;%d;%d\n'\
                       % (center_x, center_y, center_x, center_y))

        # taxon info
        from FilesIO import set_color_num
        from matplotlib import colors
        rgb_colors = [list(map(lambda x: round(x*255), c)) for c in set_color_num(len(pops))]
        colors = [i[0]+i[1]*256+i[2]*256*256 for i in rgb_colors]
        fdi_file.write('NUMBER_OF_TAXA;%d\n' % (sequence_df.index.size+len(medians))) # including median nodes
        plot_num = 0
        for name in sequence_df.index:
            if name in sample_sp_code.keys():
                size = hap_size[plot_num]
                freq = hap_freq[plot_num]
                pos_x, pos_y = list(map(round, pos[plot_num]))
                active = 'TRUE'
                plot_num += 1
            else:
                size = 1
                pos_x, pos_y = 0, 0
                active = 'FALSE'
                freq = None
            taxon_info = 'TAXON_NAME;%s;'\
                         'TAXON_NAME_HEIGHT;6;'\
                         'TAXON_NAME_COLOR;0;'\
                         'TAXON_FREQUENCY;%d;'\
                         'TAXON_ORIG_FREQUENCY;1;'\
                         'TAXON_GEOGRAPHY;;'\
                         'TAXON_PHENOTYPE;;'\
                         'TAXON_LINEAGE;;'\
                         'TAXON_GROUP1;Group1;;'\
                         'TAXON_GROUP2;Group2;;'\
                         'TAXON_GROUP3;Group3;;'\
                         'TAXON_X;%d;'\
                         'TAXON_Y;%d;'\
                         % (name, size, pos_x, pos_y)
            if freq:
                pie_info = ''
                colored_pie_num = 0
                for pie_num, pie_freq in enumerate(freq):
                    if pie_freq == 0:
                        continue
                    else:
                        colored_pie_num += 1
                        pie_color_code = colors[pie_num]
                        pie_info += 'TAXON_COLOR_PIE%d;%d;'\
                                    'TAXON_PIE_FREQUENCY%d;%d;'\
                                    'TAXON_STYLE_PIE%d;SOLID;'\
                                    % (colored_pie_num, pie_color_code, colored_pie_num, pie_freq, colored_pie_num)
            else:
                pie_info = 'TAXON_COLOR_PIE1;65535;'\
                           'TAXON_PIE_FREQUENCY1;1;'\
                           'TAXON_STYLE_PIE1;SOLID;'\

            line_info = 'TAXON_LINE_WIDTH;1;'\
                        'TAXON_LINE_COLOR;0;'\
                        'TAXON_LINE_STYLE;SOLID;'\
                        'TAXON_ACTIVE;%s\n'\
                        % (active)

            fdi_file.write(taxon_info + pie_info + line_info)

        # median info
        for n, m in enumerate(medians):
            pos_x, pos_y = list(map(round, pos[m]))
            median_info = 'TAXON_NAME;mv%d;'\
                          'TAXON_NAME_HEIGHT;6;'\
                          'TAXON_NAME_COLOR;0;'\
                          'TAXON_FREQUENCY;0;'\
                          'TAXON_ORIG_FREQUENCY;0;'\
                          'TAXON_GEOGRAPHY;;'\
                          'TAXON_PHENOTYPE;;'\
                          'TAXON_LINEAGE;;'\
                          'TAXON_GROUP1;Group1;;'\
                          'TAXON_GROUP2;Group2;;'\
                          'TAXON_GROUP3;Group3;;'\
                          'TAXON_X;%d;'\
                          'TAXON_Y;%d;'\
                          'TAXON_COLOR_PIE1;255;'\
                          'TAXON_PIE_FREQUENCY1;0;'\
                          'TAXON_STYLE_PIE1;SOLID;'\
                          'TAXON_LINE_WIDTH;1;'\
                          'TAXON_LINE_COLOR;0;'\
                          'TAXON_LINE_STYLE;'\
                          'SOLID;TAXON_ACTIVE;TRUE\n'\
                          % (n, pos_x, pos_y)
            fdi_file.write(median_info)

        fdi_file.write('\n')
        fdi_file.write('NUMBER_OF_LINKS;%d\n' % len(network_data))
        for n, link in enumerate(network_data):
            if n < torso_num:
                in_torso = 'TRUE'
            else:
                in_torso = 'FALSE'
            link_names = [list(sample_sp_code.keys())[name_code] if name_code < len(sample_sp_code) else 'mv'+str(medians.index(name_code)) for name_code in link[:2]]
            link_info = 'LINK_TAXON1;%s;'\
                        'LINK_TAXON2;%s;'\
                        'LINK_IN_TORSO;%s;'\
                        'LINK_TORSO_PART;1;'\
                        'NUMBER_OF_MUTATIONS;%d;'\
                        % (link_names[0], link_names[1], in_torso ,link[2])
            mutation_data = self._get_mutations(link, sp_sequence_df, STR_code_pos)
            mutation_info = ''
            for mutation_n, mutation in enumerate(mutation_data):
                mutation_info += 'MUTATION%d_NAME;%s;'\
                                 'MUTATION%d_NAME_HEIGHT;6;'\
                                 'MUTATION%d_NAME_COLOR;255;'\
                                 'LINK_LINE_WIDTH;1;'\
                                 'LINK_LINE_STYLE;SOLID;'\
                                 'LINK_LINE_COLOR;8421504;'\
                                 'LINK_ACTIVE;TRUE'\
                                 % (mutation_n+1, mutation, mutation_n+1, mutation_n+1)
            fdi_file.write(link_info+mutation_info+'\n')

        # other info
        # median nodes not included
        fdi_file.write('\n')
        other_info = 'NUMBER_OF_HELP_POINT_LINKS;0\n'\
                     '\n'\
                     'LONG_EDGES;0\n'\
                     'TOTAL_TAXA;%d\n'\
                     'TOTAL_HAPLOTYPES;%d\n'\
                     % (sequence_df.index.size, len(sample_sp_code))
        fdi_file.write(other_info)

        # equal info
        equal_info = ''
        equal_num = 0
        for sample, hap in sample_sp_code.items():
            equal_samples = sequence_df[sequence_df['Haplotype']==hap].index.tolist()
            for equal in equal_samples:
                if equal != sample:
                    equal_num += 1
                    equal_info += 'EQUIVALENT_TAXA;%s;%s\n' % (sample, equal)
        fdi_file.write('NUMBER_EQUIVALENT_TAXA;%d\n' % equal_num)
        fdi_file.write(equal_info)


        # end info
        end_info = 'COLOR_SCHEME;SCHEME_ACTIVE;FALSE;SCHEME_TYPE;Phenotype;\n'\
                   'COLOR_SCHEME;SCHEME_ACTIVE;FALSE;SCHEME_TYPE;Geography;\n'\
                   'COLOR_SCHEME;SCHEME_ACTIVE;FALSE;SCHEME_TYPE;Lineage;\n'\
                   'COLOR_SCHEME;SCHEME_ACTIVE;FALSE;SCHEME_TYPE;Group1;\n'\
                   'COLOR_SCHEME;SCHEME_ACTIVE;FALSE;SCHEME_TYPE;Group2;\n'\
                   'COLOR_SCHEME;SCHEME_ACTIVE;FALSE;SCHEME_TYPE;Group3;\n'
        fdi_file.write(end_info)

        fdi_file.close()

    # summarize parameters for figure
    def _get_figure_parameters(self, sequence_df, population_data, network_data, alt_data, node_num, original_network_size, alt_type, euqal_size, label):

        haplotypes = sorted(list(set(sequence_df['Haplotype'])))
        # calculate frequency of each haplotypes
        if euqal_size:
            hap_size = [1]*len(haplotypes)
        else:
            hap_size = [sequence_df['Haplotype'].tolist().count(hap) for hap in haplotypes]
            hap_size_per = [size / sequence_df.index.size for size in hap_size ]


        if isinstance(population_data, pd.DataFrame):
            pops = sorted(list(set(population_data[label])))
            # calculate haplotype frequency in each population
            hap_pop_data = pd.DataFrame(index=range(len(haplotypes)), columns=pops)
            for pop in pops:
                pop_inds = population_data[population_data[label]==pop]['SampleID'].tolist()
                pop_data = sequence_df.loc[pop_inds]
                hap_pop_data[pop] = [pop_data['Haplotype'].tolist().count(hap) for hap in haplotypes]
            hap_freq = [hap_pop_data.loc[i].tolist() for i in hap_pop_data.index]
            hap_freq_per = [(hap_pop_data.loc[i]/sum(hap_pop_data.loc[i])).round(6).tolist() for i in hap_pop_data.index]
        else:
            hap_freq_per = None
            pops = None

        # set length of each node
        length_data = np.empty(shape=(0, 3))
        num = np.max(network_data[:, :2])+1
        to_be_deleted_row = []
        for row_index, i in enumerate(network_data):
            step = i[2]
            if step > 1:
                start = i[0]
                end = i[1]
                for i in range(step):
                    if i == 0:
                        length = [start, num, 1]
                    elif i == step-1:
                        length = [num, end, 1]
                    else:
                        length = [num, num+1, 1]
                        num += 1
                    length_data = np.concatenate((length_data, np.array([length])))
                num += 1
                to_be_deleted_row.append(row_index)
                original_network_size -= 1
        network_data = np.delete(network_data, to_be_deleted_row, axis=0)
        length_data = length_data.astype(int)

        # set figure size and network class
        import warnings
        warnings.filterwarnings('ignore')
        font_files = font_manager.findSystemFonts(fontpaths=[os.path.split(os.path.realpath(__file__))[0] + '/sans-serif'])
        font_list = font_manager.createFontList(font_files)
        font_manager.fontManager.ttflist.extend(font_list)
        plt.rcParams['font.family']= 'sans-serif'
        plt.rcParams['font.sans-serif'] = ['Arial'] # font
        plt.rcParams['legend.handlelength'] = 1
        plt.rcParams['legend.handleheight'] = 1
        fig_len = len(length_data)/3
        fig, ax = plt.subplots(figsize=(fig_len, fig_len))
        trans = ax.transData.transform
        trans2 = fig.transFigure.inverted().transform
        DG = nx.DiGraph()

        # construct network by adding nodes and edges
        style_list = []
        for i, j in enumerate(network_data):
            if i < original_network_size:
                DG.add_edge(j[0], j[1])
                style_list.append('solid')
            else:
                DG.add_edge(j[0], j[1])
                if alt_type != 'dashed':
                    style_list.append('dashed')
                else:
                    style_list.append('solid')
        # construct network by adding node of length
        for i in length_data:
            DG.add_edge(i[0], i[1])
            style_list.append('solid')

        try:
            import pygraphviz as pg
            pos = nx.nx_agraph.graphviz_layout(DG, prog='neato')
        except ImportError:
            pre_pos = nx.spring_layout(DG)
            pos = {}
            for i in pre_pos.keys():
                pos[i] = trans(pre_pos[i])

        return length_data, DG, pos, hap_size, hap_freq, hap_size_per, hap_freq_per, style_list, pops, trans, trans2, fig_len

    # plot a tree and replace nodes to pies
    def plot_network(self, length_data, DG, pos, node_num, median_size_rate, hap_size_per, hap_freq_per, style_list, pops, trans, trans2, fig_len, label):


        import warnings
        warnings.filterwarnings('ignore')

        # set figure parameters
        # mdian haplotype size
        base_size = 5*20
        median_size = base_size * median_size_rate
        # all node size
        hap_num = len(hap_size_per)
        max_size = max(hap_size_per)
        size_rate = base_size/max_size
        size = [hap_size_per[i]*size_rate if i < hap_num else median_size for i in DG.nodes][:node_num]
        size = size + [0.2] * length_data.size

        try:
            # color and line type
            if hap_freq_per:
                from matplotlib.patches import Patch
                from FilesIO import set_color_num
                nx.draw(DG, pos=pos, arrows=False, node_color='black', node_size=min(size), style=style_list)
                # set colors
                colors = set_color_num(len(pops))
                legend_elements = []
                for color, pop in zip(colors , pops):
                    legend_elements.append(Patch(facecolor=color, edgecolor='black', label=pop))
                plt.legend(handles=legend_elements, frameon=False, bbox_to_anchor=(-0.04, 0.5), prop={'size': fig_len})
                pie_size = 0.025
                move = pie_size/2.0
                for i in DG:
                    if i < hap_num:
                        xx, yy = trans(pos[i])
                        xa, ya = trans2((xx,yy))
                        pie_ax = plt.axes([xa-move, ya-move, pie_size, pie_size])
                        pie_ax.set_aspect('equal')
                        fracs = hap_freq_per[i]
                        if 1 in fracs:
                            one_colors= [colors[fracs.index(1)]]
                            one_fracs = [1]
                            pie_ax.pie(one_fracs, colors=one_colors, radius=(1/(max_size/hap_size_per[i]))**0.5, wedgeprops={'edgecolor':'black','linewidth': 0.5, 'linestyle': 'solid'})
                        else:
                            pie_ax.pie(fracs, colors=colors, radius=(1/(max_size/hap_size_per[i]))**0.5, wedgeprops={'edgecolor':'black','linewidth': 0.5, 'linestyle': 'solid'})
            else:
                colors = ['yellow' if i < hap_num else 'black' for i in DG.nodes] # color
                nx.draw(DG, pos=pos, arrows=False, node_color=colors, node_size=size, style=style_list, linewidths=1)
                node_ax = plt.gca() # to get the current axis
                node_ax.collections[0].set_edgecolor('black')
                node_ax.collections[0].set_linewidth(0.5)

            # output
            if label:
                plt.savefig(self.path+'.net.fig.%s.pdf' % label, dpi=100, bbox_inches='tight')
            else:
                plt.savefig(self.path+'.net.fig.pdf', dpi=100, bbox_inches='tight')
            plt.close()
        except AttributeError:
            print('[Y-LineageTracker] Cannot generate figure since matplotlib.cbook version bugs')

    # main function of network analysis
    def calculate_and_output(self, haplotype_data, population_data, arguments, label):

        # convert data to 0-1 code and calculate distance matrix
        sequence_df, STR_code_pos = self._convert_allel_to_binary(haplotype_data, arguments.filter)
        distance_matrix, sp_sequence_df, sample_sp_code = self._calculate_distance_matrix(sequence_df)
        # calculate network
        self.logger.info('[Y-LineageTracker] Calculating median haplotypes')
        network_data, alt_data, sp_sequence_df2 = self._calculate_network(distance_matrix, sp_sequence_df, arguments.treetype, arguments.search)
        plot_network_data, node_num = self._get_network_for_plot(network_data, alt_data, sp_sequence_df.index.size, arguments.alt)

        # visulization
        self.logger.info('[Y-LineageTracker] Start to plot network')

        if label != 'Group':
            length_data, DG, pos, hap_size, hap_freq, hap_size_per, hap_freq_per, style_list, pops, trans, trans2, fig_len = self._get_figure_parameters(sequence_df, population_data, plot_network_data, alt_data, node_num, network_data.size, arguments.alt, arguments.euqalsize, label)
            if arguments.fdi:
                self._get_fdi_file(plot_network_data, plot_network_data.size, node_num, sequence_df, sp_sequence_df2, sample_sp_code, pos, hap_size, hap_freq, STR_code_pos, pops)
            self.plot_network(length_data, DG, pos, node_num, arguments.mediansize, hap_size_per, hap_freq_per, style_list, pops, trans, trans2, fig_len, label)
        else:
            label = ['Population', 'Group']
            for i in label:
                length_data, DG, pos, hap_size, hap_freq, hap_size_per, hap_freq_per, style_list, pops, trans, trans2, fig_len = self._get_figure_parameters(sequence_df, population_data, plot_network_data, alt_data, node_num, network_data.size, arguments.alt, arguments.euqalsize, label)
                if arguments.fdi:
                    self._get_fdi_file(plot_network_data, plot_network_data.size, node_num, sequence_df, sp_sequence_df2, sample_sp_code, pos, hap_size, hap_freq, STR_code_pos, pops)
                self.plot_network(length_data, DG, pos, node_num, arguments.mediansize, hap_size_per, hap_freq_per, style_list, pops, trans, trans2, fig_len, label)

        self.logger.info('[Y-LineageTracker] Network analysis finished')



def main():

    start = time.perf_counter()
    arguments = net_parser()

    from FilesIO import check_seq_matrix_input, check_hap_input, set_haplotype_index, check_population, check_overlap, check_pop_header, time_count
    # check input files
    input, format = check_seq_matrix_input(arguments.seq, arguments.format, arguments.matrix)

    # check haplotype data
    haplotype_data = check_hap_input(input, 'haplotype', format)

    # check filter argument
    check_misisng_rate(arguments.filter)

    from FilesIO import get_out_path
    # get path of output files
    path = get_out_path(arguments.output, 'Net')
    log_file = path + '.NetLog.log'

    # set log file
    set_log(arguments, log_file)

    # check if there is a population file
    if arguments.population:
        population_data = check_population(arguments.population)
        population_data = check_pop_header(population_data)
        haplotype_data, population_data = check_overlap(haplotype_data, population_data)
        label = 'Population'
        if 'Group' in population_data.columns:
            label = 'Group'
    else:
        population_data = None
        label = ''

    # set index of haplotype data
    haplotype_data = set_haplotype_index(haplotype_data, format)

    # calculate network and visualize result
    net = HaploNetwork(arguments.gap, path)
    net.calculate_and_output(haplotype_data, population_data, arguments, label)

    time_count(start)


if __name__ == '__main__':
    main()
