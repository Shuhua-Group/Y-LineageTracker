import os
import re
import sys
import time
import logging
import argparse
import pandas as pd
import numpy as np


'''
Description:
This module is used for estimating TMRCA of haplotypes.
Network information file or sequence alignment file is supported to estimate TMRCA.
The program will calculate rho statistics from network information file, or calculate ASD statistics from the sequence alignment file.
Then estimate TMRCA according to statistics.
'''


def tmrca_parser():

    parser = argparse.ArgumentParser('tmrca', description='(c) Y-LineageTracker: TMRCA estimation')
    # function used for TMRCA estimation
    parser.add_argument('tmrca',
                        help='Estimate TMRCA of Y haplotypes.')
    # required, network information file or matrix file
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--net',
                        required=False,
                        type=str,
                        action='store',
                        help='net: Network information file.')
    group.add_argument('--matrix',
                        required=False,
                        type=str,
                        action='store',
                        help='matrix: Y-STR haplotype data in matrix format')
    # required, information of ancestor and descendants.
    parser.add_argument('--anc-des',
                        required=True,
                        type=str,
                        dest='ad',
                        action='append',
                        help='anc-des: Information of ancestor and descendants. \
                              The program will estimate the time based on mutations from the ancestor to each descendant.')
    # required, population file
    parser.add_argument('-p', '--population',
                        required=False,
                        type=str,
                        action='store',
                        help='population: A file containing sample ID and population information of each individual.')
    # optional, mutation rate of Y-STR
    parser.add_argument('--mut-rate',
                        required=False,
                        type=float,
                        dest='rate',
                        action='store',
                        default=6e-3,
                        help='mut-rate: The general mutation rate across Y-STR haplotype.')
    # optional, generation time
    parser.add_argument('--generation',
                        required=False,
                        type=int,
                        action='store',
                        default=30,
                        help='generation: The generation time used for time estimation.')
    # optional, the prefix of output
    parser.add_argument('-o', '--output',
                        required=False,
                        type=str,
                        action='store',
                        help='output: The prefix of output files.')

    args = parser.parse_args()

    return args


# print program information and write to log file
def set_log(log_file, args_log):

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

    log_info = ['[Y-LineageTracker] [TMRCA]',
                '[Y-LineageTracker] Run Date: ' + time.asctime(time.localtime(time.time()))]

    if args_log.net:
        log_info.append('[Y-LineageTracker] Network Info File: %s' % args_log.net)
        stat = 'Rho'
    elif args_log.matrix:
        log_info.append('[Y-LineageTracker] Matrix File: %s' % args_log.matrix)
        stat = 'ASD'

    log_info.extend(['[Y-LineageTracker] Ancestor and Descendants Information: %s' % args_log.ad,
                     '[Y-LineageTracker] Statistics: %s' % stat,
                     '[Y-LineageTracker] Mutation Rate: %s' % args_log.rate,
                     '[Y-LineageTracker] Generation: %s' % args_log.generation])

    print('\n')
    for i in log_info:
        logger.info(i)


class CommonAncestor(object):
    '''
    This class is used to estimate TMRCA of haplotypes
    The program can calculate:
    1. rho statistics from network information file
    2. ASD statistics from matrix file
    There is only one output, that is the TMRCA of haplotypes
    Note that the time estimated by ASD is larger than that of rho
    '''

    def __init__(self):

        from GetConfig import getConfig
        config = getConfig()
        self.logger = logging.getLogger()
        self.skip_letters = config.get('Statistics', 'SkipLetters').split(',')

    def _parse_hap_data(self, matrix_file, population_file):

        from FilesIO import check_hap_input, check_population, check_pop_header, check_overlap, set_haplotype_index
        hap_data = check_hap_input(matrix_file, 'haplotype', 'matrix')

        population_data = check_population(population_file)
        population_data = check_pop_header(population_data)

        hap_data, population_data = check_overlap(hap_data, population_data)
        hap_data = set_haplotype_index(hap_data, 'matrix')
        population_data = population_data.set_index(['SampleID'], drop=True)

        return hap_data, population_data

    def _parse_net_data(self, net_file):

        net = open(net_file, 'r').readlines()

        for i, j in enumerate(net):
            if j.startswith('#HapLength'):
                length = int(j.strip().split(':')[-1])
            elif j.startswith('#HapSamples'):
                sample_line = i
            elif j.startswith('#HapLinks'):
                link_line = i

        sample_info = net[sample_line+1:link_line-1]
        sample_num = {}
        for i in sample_info:
            line = i.strip()
            num = int(line.split('\t')[0])
            samples = line.split('\t')[-1].split(',')
            for j in samples:
                sample_num[j] = num

        link_info = net[link_line+1:]
        link_data = np.array([i.strip().split('\t') for i in link_info])
        link_data = link_data.astype(int)

        net_data = [length, sample_num, link_data]

        return net_data

    def _get_ancestor_and_descendants(self, anc_des_list):

        if len(anc_des_list) == 1:
            anc_des = anc_des_list[0]
            if os.path.isfile(anc_des):
                pre_ad_list = open(anc_des).readlines()
            else:
                pre_ad_list = [anc_des]
        else:
            pre_ad_list = anc_des_list

        ad_list = []
        for i in pre_ad_list:
            if re.match(r'^\w+:(\w+,)+\w+$', i):
                ad_list.append(i)
            else:
                self.logger.error('[Y-LineageTracker] %s format error' % i)
                sys.exit()

        return ad_list

    def _calculate_tmrca(self, data, rate, generation, ad_list, stat, path, population_data):

        if stat == 'rho':
            L = data[0]
            sample_num = data[1]
            link_data = data[2]
        elif stat == 'asd':
            L = data.columns.size # haplotype length
        mu = rate # mutation rate
        g = generation # generation

        num = 0
        out_df = pd.DataFrame(columns=['Ancestor', 'Descendants', stat, 'TMRCA'])

        for ad in ad_list:
            a = ad.split(':')[0]
            d = ad.split(':')[1].split(',')

            if stat == 'rho':
                ancestor_node = sample_num[a]
                descendant_nodes = [sample_num[i] for i in d]
                res = self._calculate_rho(ancestor_node, descendant_nodes, link_data)
            elif stat == 'asd':
                pop_inds = {}
                for pop in sorted(set(population_data['Population'])):
                    inds = population_data[population_data['Population']==pop].index.tolist()
                    pop_data = data.loc[inds].to_numpy().T
                    pop_inds[pop] = pop_data

                ancestor_hap = data.loc[a, :].tolist()
                ancestor_pop = population_data.at[a, 'Population']
                ancestor_data = pop_inds[ancestor_pop]
                ancestor_freq = [list(j).count(i) / len(list(filter(lambda x: x not in self.skip_letters, j))) for i, j in zip(ancestor_hap, ancestor_data)]

                descendant_haps = [data.loc[i, :].tolist() for i in d]
                descendant_freqs = []
                for d_ind, d_hap in zip(d, descendant_haps):
                    descendant_pop = population_data.at[d_ind, 'Population']
                    descendant_data = pop_inds[descendant_pop]
                    descendant_freq = [list(j).count(i) / len(list(filter(lambda x: x not in self.skip_letters, j))) for i, j in zip(d_hap, descendant_data)]
                    descendant_freqs.append(descendant_freq)

                res = self._calculate_asd(ancestor_hap, descendant_haps, ancestor_freq, descendant_freqs)
            time = round((res / (mu*L)) * g, 6)
            line = [a, ','.join(d), res, time]
            out_df.loc[num] = line
            num += 1

        out_df.to_csv(path + '.%s.txt' % stat, sep='\t', index=False)

    # the average number of sites differing between a set of haplotypes and their specified common ancestor
    def _calculate_asd(self, ancestor_hap, descendant_haps, ancestor_freq, descendant_freqs):

        diffs = []
        for d_hap, d_freq in zip(descendant_haps, descendant_freqs):
            diff_num = np.sum([(int(i)-int(j))**2*k*l for i, j, k, l in zip(ancestor_hap, d_hap, ancestor_freq, d_freq) if i not in self.skip_letters and j not in self.skip_letters])
            diffs.append(diff_num)

        res = round(np.sum(diffs), 2)

        return res

    def _calculate_rho(self, ancestor_node, descendant_nodes, link_data):

        end = ancestor_node
        diff_num = []
        for d_node in descendant_nodes:
            start = d_node
            diff = self._calculate_net_rho(start, end, link_data)
            diff_num.append(diff)

        res = round(np.mean(diff_num), 2)

        return res

    def _calculate_net_rho(self, start, end, link_data):

        vs = [start]
        num = 0
        while True:
            if num == 0:
                vs2 = [start]
                steps = [0]
            else:
                vs2 = vs3
                steps = steps2
            vs3 = []
            steps2 = []
            for i, j in zip(vs2, steps):
                v = link_data[np.where(link_data[:, :2]==i)[0]]
                nodes = sorted(list(set(v[:,:2].reshape(-1))))
                nodes2 = list(filter(lambda x: not x in vs, nodes))
                stp = (v[:,2]+j).tolist()
                stp = [stp[np.where(v==n)[0][0]] for n in nodes2]
                vs3.extend(nodes2)
                steps2.extend(stp)
                vs.extend(vs3)
            num += 1
            if end in vs3:
                diff_num = min(steps2)
                break

        return diff_num

    def run_tmrca(self, arguments, path):

        if arguments.net:
            data = self._parse_net_data(arguments.net)
            stat = 'rho'
            population_data = None
        elif arguments.matrix:
            if arguments.population:
                data, population_data = self._parse_hap_data(arguments.matrix, arguments.population)
                stat = 'asd'
            else:
                self.logger.error('[Y-LineageTracker] Population file is required for matrix file to perform ASD estimation')
                sys.exit()

        self.logger.info('[Y-LineageTracker] Reading input data')
        ad_list = self._get_ancestor_and_descendants(arguments.ad)
        self.logger.info('[Y-LineageTracker] Get ancestor and descendants')
        self.logger.info('[Y-LineageTracker] Calculating TMRCA')
        self._calculate_tmrca(data, arguments.rate, arguments.generation, ad_list, stat, path, population_data)
        self.logger.info('[Y-LineageTracker] Finished')


def main():

    from FilesIO import get_out_path, time_count
    start = time.perf_counter()
    arguments = tmrca_parser()
    path = get_out_path(arguments.output, 'TMRCA')
    log_file = path + '.TMRCALog.log'

    set_log(log_file, arguments)

    estimate = CommonAncestor()
    estimate.run_tmrca(arguments, path)

    time_count(start)

if __name__ == '__main__':
    main()
