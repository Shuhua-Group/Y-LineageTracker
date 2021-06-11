import sys
import time
import random
import logging
import argparse
import pandas as pd
import numpy as np

from itertools import chain, repeat
from collections import Counter
from GetConfig import getConfig


config = getConfig()


'''
Description:
This module is used for statistical analysis.
Y haplotype data is required for analysis.
It can calculate commonly used statistics in population genetics.
HD: haplotype diversity
MPD: mean pairwise differences
GD: genetic distance
AMOVA: analysis of molecular variance
predict: predict haplogroup from Y-STR haplotype
'''


def stat_parser():

    parser = argparse.ArgumentParser('stat', description='(c) Y-LineageTracker: Statistical analysis')
    # function used for statistical analysis
    parser.add_argument('stat',
                        help='Perform statistical analysis from Y haplotype data.')
    # required, input file, alignment sequence or matrix
    input = parser.add_mutually_exclusive_group(required=True)
    input.add_argument('--seq',
                        type=str,
                        action='store',
                        help='seq: Haplotype data in sequence alignment format.')
    input.add_argument('--matrix',
                        type=str,
                        action='store',
                        help='matrix: Y-STR haplotype data in matrix format.')
    # optional, format of sequence alignment
    parser.add_argument('--seq-format',
                        required=False,
                        type=str,
                        dest='format',
                        action='store',
                        choices=['fasta', 'phylip', 'nexus', 'meg', 'vcf'],
                        help='seq-format: The format of sequence file. \
                              This option is only required for the sequence alignment file.')
    # population file, required for calculating statistics, not required for prediction
    parser.add_argument('-p', '--population',
                        required=False,
                        type=str,
                        action='store',
                        help='population: A file containing sample ID and population information of each individual. \
                              The population file is required in the analysis of calculating statistics, but not in the haplogroup prediction')
    # optional, calculate genetic distance
    parser.add_argument('--gd',
                        required=False,
                        type=str,
                        action='store',
                        choices=['fst', 'rst', 'Ngst', 'NCgst', 'Hgst', 'ibs'],
                        help='gd: Calculate genetic distance based on the selected statistics.')
    # optional, perform analysis of molecular variance
    parser.add_argument('--amova',
                        required=False,
                        type=str,
                        action='store',
                        nargs='?',
                        choices=['fst', 'rst'],
                        help='amova: Get the AMOVA table and calculate Fst or Rst of pairwise populations based on AMOVA and give the p values. \
                              The default is fst.')
    # optional, calculate haplotype diversity
    parser.add_argument('--hd',
                        required=False,
                        action='store_true',
                        help='hd: Calculate haplotype diversity of each population.')
    # optional, calculate mean pairwise differences
    parser.add_argument('--mpd',
                        required=False,
                        action='store_true',
                        help='mpd: calculate mean pairwise differences within and between populations.')
    # optional, perform prediction anlaysis from Y-STR haplotype
    parser.add_argument('--predict',
                        required=False,
                        action='store_true',
                        help='predict: predict possible NRY haplogroups and give the probability of each haplogroup from Y-STR haplotype data by Bayesian approach. \
                              This analysis only supports Y-STR haplotype data.')
    # optional, the prefix of output
    parser.add_argument('-o', '--output',
                        required=False,
                        type=str,
                        action='store',
                        help='output: The prefix of output files.')

    args = parser.parse_args()

    return args


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

    log_info = ['[Y-LineageTracker] [Stat]',
                '[Y-LineageTracker] Run Date: ' + time.asctime(time.localtime(time.time())),
                '[Y-LineageTracker] Haplotype File: %s' % args_log.seq,
                '[Y-LineageTracker] Haplotype File Format: %s' % args_log.format]

    if args_log.mpd:
        log_info.append('[Y-LineageTracker] Calculate mean pairwise differences')
    if args_log.hd:
        log_info.append('[Y-LineageTracker] Calculate haplotype diversity')
    if args_log.amova:
        log_info.append('[Y-LineageTracker] Calculate %s based on AMOVA' % args_log.amova)
    if args_log.predict:
        log_info.append('[Y-LineageTracker] Predict haplogroup from Y-STR data')

    print('\n')
    for i in log_info:
        logger.info(i)


# check input command
def check_commands(arguments):

    if arguments.mpd or arguments.hd or arguments.amova or arguments.gd:
        if not arguments.population:
            print('[Y-LineageTracker] Population file required for specified analysis')
            sys.exit()
    else:
        if not arguments.predict:
            print('[Y-LineageTracker] No command for analysis')
            sys.exit()


# check amova command
def check_amova(amova, arguments):

    if amova == None:
        if '--amova' in arguments:
            amova = 'fst'
    else:
        amova = amova

    return amova


class StatisticsAnalysis(object):
    '''
    This class include commonly used statistical method for Y-chromosome analysis:
    1. HD: haplotype diversity
    2. MPD: mean pairwise differences
    3. GD: genetic distance
    4. AMOVA: analysis of molecular analysis
    for analysis 2-4, p value will calculated by permutation
    '''

    def __init__(self, label, haplotype_data, population_data):

        self.logger = logging.getLogger()
        self.label = label
        if label == 'Population':
            self.label_num = 1
        else:
            self.label_num = 2
        self.skip_letters = config.get('Statistics', 'SkipLetters').split(',')
        self.haplotype_array = np.array(haplotype_data)
        self.haplotype_array_with_inds = haplotype_data.reset_index().to_numpy()
        self.population_array = np.array(population_data)
        self.header = haplotype_data.columns.tolist()
        self.pops = sorted(list(set(population_data[label])))

    # calculate distance between two sequence
    def _pairwise_distance(self, seq1, seq2, type):

        if type == 'num_of_diff_alleles':
            distance = sum(l1 != l2 for l1, l2 in zip(seq1, seq2) if l1 not in self.skip_letters and l2 not in self.skip_letters)
        elif type == 'sum_of_squared_diff':
            distance = sum((int(l1)-int(l2))**2 for l1, l2 in zip(seq1, seq2) if l1 not in self.skip_letters and l2 not in self.skip_letters)

        return distance

    # get haplotype data of a population
    def _get_population_data(self, pop, population_data=None):

        # string type: one population, list type: two populations or more
        if isinstance(population_data, np.ndarray):
            population_data = population_data
        else:
            population_data = self.population_array
        if isinstance(pop, str):
            pop_inds = population_data[population_data[:, self.label_num]==pop][:, 0]
        elif isinstance(pop, list):
            pop_inds = population_data[list(map(lambda x: x in pop, population_data[:, self.label_num]))][:, 0]
        pop_haplotype_data = self.haplotype_array[list(map(lambda x: x in pop_inds, self.haplotype_array_with_inds[:, 0]))]

        return pop_haplotype_data

    # calculate haplotype frequency of a population
    def _calculate_haplotype_freq(self, pop_haplotype_data):

        pop_all_haps = [tuple(i) for i in pop_haplotype_data] # get all haplotypes
        pop_haps = sorted(list(set(pop_all_haps))) # get haplotypes without duplication
        hap_freq = [pop_all_haps.count(i) / len(pop_all_haps) for i in pop_haps] # calculate freq of each haplotype for the population

        return pop_haps, hap_freq

    # calculate allele frequency of a population
    def _calculate_allele_freq(self, pops, population_data, pops_haplotype_data):

        # get allele frequency, sample size of each population
        pops_count = []
        sample_sizes = []
        max_len = max([len(set(pops_haplotype_data[:, i])) for i in range(len(self.header))])
        for pop in pops:
            pop_haplotype_data = self._get_population_data(pop, population_data)
            single_pop_count = []
            sample_sizes.append(len(pop_haplotype_data))
            for allele in range(len(self.header)):
                all_alleles = sorted(list(set(pops_haplotype_data[:, allele])))
                allele_count = [list(pop_haplotype_data[:, allele]).count(i) for i in all_alleles]
                if len(allele_count) < max_len:
                    diff = max_len - len(allele_count)
                    allele_count += [np.nan] * diff
                single_pop_count.append(allele_count)
            pops_count.append(single_pop_count)
        pops_count = np.array(pops_count)
        sample_sizes = np.array(sample_sizes)

        return pops_count, sample_sizes

    # calculate haplotype diversity of a population
    def _haplotype_diversity(self, pop, type, allele=None):

        pop_haplotype_data = self._get_population_data(pop)
        n = len(pop_haplotype_data)
        if type == 'single':
            hap_freq = (np.unique(pop_haplotype_data[:, allele], return_counts=True)[1]/len(pop_haplotype_data[:, allele]))
        elif type == 'haplotype':
            pop_haps, hap_freq = self. _calculate_haplotype_freq(pop_haplotype_data)

        if n == 1:
            return 'NA', 'NA'

        # calculate single Y-STR and Y-STR haplotype diversity
        sub = 0
        for i in hap_freq:
            sub += i**2
        diversity = (n/(n-1)) * (1-sub)

        # calculate var and sd
        sub1 = 0
        sub2 = 0
        for i in hap_freq:
            sub1 += i**3
            sub2 += i**2

        v1 = 2.0 / (n*(n-1))
        v2 = 2.0 * (n-2) * (sub1-sub2**2)

        var = v1 * (v2+sub2-sub2**2)
        sd = var**0.5

        return diversity, sd

    # get haplotype diversity of all populations
    def _calculate_haplotype_diversity(self):

        hd_df = pd.DataFrame(index=self.pops, columns=['Haplotype']+self.header)
        for pop in self.pops:
            for n, allele in enumerate(self.header):
                diversity, sd = self._haplotype_diversity(pop, 'single', n)
                hd_df.at[pop, allele] = str(round(diversity, 5))+'±'+str(round(sd, 5))
            diversity, sd = self._haplotype_diversity(pop, 'haplotype')
            if diversity == 'NA':
                hd_df.at[pop, 'Haplotype'] = 'NA'
            else:
                hd_df.at[pop, 'Haplotype'] = str(round(diversity, 5))+'±'+str(round(sd, 5))

        return hd_df

    # calculate mean pairwise differences of between two populations
    def _between_mean_pairwise_differences(self, pop1, pop2):

        pop1_data = self._get_population_data(pop1)
        pop2_data = self._get_population_data(pop2)
        n = len(pop1_data) + len(pop2_data)
        pop1_haps, hap1_freq = self._calculate_haplotype_freq(pop1_data)
        pop2_haps, hap2_freq = self._calculate_haplotype_freq(pop2_data)

        pi = 0
        for i in range(len(pop1_haps)):
            Sub = 0
            for j in range(len(pop2_haps)):
                hap1 = pop1_haps[i]
                hap2 = pop2_haps[j]
                d = self._pairwise_distance(hap1, hap2, 'num_of_diff_alleles')
                p_i = hap1_freq[i]
                p_j = hap2_freq[j]
                Sub += p_i*p_j*d
            pi += Sub

        if n > 6:
            # calculate var ans sd
            sub1 = 3 * n * (n+1)
            sub2 = 2 * (n**2+n+3)
            sub3 = 11 * (n**2-7*n+6)

            var = float(sub1*pi+sub2*pi**2) / sub3
            sd = var**0.5
        else:
            sd = 'NA'

        return pi, sd

    # calculate mean pairwise difference within a population
    def _within_mean_pairwise_differences(self, pop):

        pop_haplotype_data = self._get_population_data(pop)
        n = len(pop_haplotype_data)
        pop_haps, hap_freq = self. _calculate_haplotype_freq(pop_haplotype_data)

        # n: sample size
        # hap_freq: frequence of haplotypes of two population
        if round(sum(hap_freq), 1) != 1.0:
            return 'NA', 'NA'
        if n <= 2:
            return 'NA', 'NA'

        pi = 0
        k = len(pop_haps)
        for i in range(k):
            Sub = 0
            for j in range(k):
                hap1 = pop_haps[i]
                hap2 = pop_haps[j]
                d = self._pairwise_distance(hap1, hap2, 'num_of_diff_alleles')
                p_i = hap_freq[i]
                p_j = hap_freq[j]
                Sub += p_i*p_j*d
            pi += Sub
        pi = (float(n)/(n-1)) * pi

        if n > 6:
            # calculate var ans sd
            sub1 = 3 * n * (n+1)
            sub2 = 2 * (n**2+n+3)
            sub3 = 11 * (n**2-7*n+6)
            var = float(sub1*pi+sub2*pi**2) / sub3
            sd = var**0.5
        else:
            sd = 'NA'

        return pi, sd

    # get mean pairwise difference of all populations
    def _calculate_mean_pairwise_differences(self):

        mpd_df = pd.DataFrame(index=self.pops, columns=self.pops)
        pop_num = 0
        for pop1 in self.pops:
            for pop2 in self.pops[pop_num:]:
                if pop1 == pop2:
                    pi, sd = self._within_mean_pairwise_differences(pop1)
                else:
                    pi, sd = self._between_mean_pairwise_differences(pop1, pop2)
                if pi == 'NA':
                    mpd_df.at[pop2, pop1] = 'NA'
                else:
                    if sd == 'NA':
                        mpd_df.at[pop2, pop1] = str(round(pi, 5))+'±'+'NA'
                    else:
                        mpd_df.at[pop2, pop1] = str(round(pi, 5))+'±'+str(round(sd, 5))
            pop_num += 1

        return mpd_df

    # claculate pairwise Fst
    # Weir, B. S. and Hill, W. G. (2002) Estimating F-statistics. Annual Review of Genetics, 36, 721–750.
    def _fst(self, pops, population_data, pops_haplotype_data):

        # get allele frequency, sample size of each population
        pops_count, sample_sizes = self._calculate_allele_freq(pops, population_data, pops_haplotype_data)

        # calculate MSP and MSG to get theta statistic
        pops_freq = [np.array(i) / j for i, j in zip(pops_count, sample_sizes)]
        mean_freq = np.sum([np.array(i) for i in pops_count]) / np.sum(sample_sizes)


        MSP = np.sum([i*((j-mean_freq)**2) for i, j in zip(sample_sizes, pops_freq)], axis=0)

        MSG1 = float(1/np.sum([i-1 for i in sample_sizes]))
        MSG2 = np.sum([np.array(i)*(1-j) for i, j in zip(pops_count, pops_freq)], axis=0)
        MSG = MSG1 * MSG2

        nc = np.sum(sample_sizes) - np.sum(sample_sizes**2)/np.sum(sample_sizes)

        theta1 = np.nansum(MSP)-np.nansum(MSG)
        theta2 = np.nansum(MSP)+(nc-1)*np.nansum(MSG)
        if theta2 == 0:
            theta = 0
        else:
            theta = (theta1 / theta2)

        return theta

    # calculate pairwise Rst
    # Slatkin, M. (1995) A measure of population subdivision based on microsatellite allele frequencies. Genetics, 139, 457–462.
    def _rst(self, pops, population_data, distance_matrix):

        n = len(distance_matrix)
        pop_num = len(pops)
        expr = 2 / (pop_num*(pop_num-1))

        level_labels = population_data[0:, self.label_num]
        labels_num = []
        for i in Counter(level_labels).values():
            labels_num += [i]*i
        level_labels_code = [sorted(list(set(level_labels))).index(i) for i in level_labels]

        # calculate SW
        SW = 0
        bool1 = np.equal.outer(level_labels_code, level_labels_code)
        for i in range(pop_num):
            bool2 = np.tile(np.array(level_labels_code) == 1, (n, 1))
            bool_SW = bool1 & bool2 & bool2.T
            SW += np.sum(distance_matrix[bool_SW]) / np.sum(bool_SW)
        SW = SW / pop_num

        # calculate SB
        SB = 0
        bool3 = np.not_equal.outer(level_labels_code, level_labels_code)
        for i in range(pop_num-1):
            for j in range(i+1, pop_num):
                bool4 = np.tile(np.array(level_labels_code) == i, (n, 1)) | np.tile(np.array(level_labels_code) == j, (n, 1))
                bool_SB = bool3 & bool4 & bool4.T
                SB += np.sum(distance_matrix[bool_SB]) / np.sum(bool_SB)
        SB = SB * expr

        n_bar = np.mean(list(Counter(level_labels).values()))
        denom = np.sum(list(Counter(level_labels).values())) - 1
        S_bar = (SW*(n_bar-1)/denom) + (SB*(n_bar*(pop_num-1))/denom)
        res = (S_bar-SW) / S_bar

        return res

    # calculate hs  standard hs and ht of Gst
    def _gst(self, pops, population_data, pops_haplotype_data):

        # get allele frequency, sample size of each population
        pops_count, sample_sizes = self._calculate_allele_freq(pops, population_data, pops_haplotype_data)

        pops_freq = np.array([np.array(i) / j for i, j in zip(pops_count, sample_sizes)])

        hs = 1 - np.nansum(np.nansum((pops_freq**2), axis=0), axis=1)/len(pops)
        ht = 1 - np.nansum((np.nansum(pops_freq**2, axis=2)/len(pops))**2, axis=0)

        return hs, ht, sample_sizes

    # calculate Nei Gst
    # Nei M (1973) Analysis of Gene Diversity in Subdivided Populations. Proc. Nat. Acad. Sci., 70, 3321-3323.
    def _Nei_gst(self, pops, population_data, pops_haplotype_data):

        hs, ht, sample_sizes = self._gst(pops, population_data, pops_haplotype_data)
        Nei_gst = 1 - (np.mean(hs)/np.mean(ht))

        return Nei_gst

    # calculate Nei and Chesser Gst
    # Nei M, Chesser RK (1983) Estimation of fixation indices and gene diversity. Annals of Human Genetics, 47, 253-259.
    def _Nei_Chesser_gst(self, pops, population_data, pops_haplotype_data, return_flag=False):

        hs, ht, sample_sizes = self._gst(pops, population_data, pops_haplotype_data)

        nm = len(pops) / np.sum(1/sample_sizes)
        Hs = (nm/(nm-1)) * hs
        Ht = ht + Hs/(nm*len(pops))

        Nei_Chesser_gst = 1 - np.mean(Hs)/np.mean(Ht)

        if return_flag:
            return Hs, Nei_Chesser_gst
        else:
            return Nei_Chesser_gst

    # calculate Hedrick Gst
    # Hedrick P (2005) A standardized genetic differentiation measure. Evolution, 59, 1633-1638.
    def _Hedrick_gst(self, pops, population_data, pops_haplotype_data):

        Hs, Nei_Chesser_gst = self._Nei_Chesser_gst(pops, population_data, pops_haplotype_data, True)

        Hedrick_gst = Nei_Chesser_gst * (len(pops)-1+np.mean(Hs))/((len(pops)-1)*(1-np.mean(Hs)))

        return Hedrick_gst

    # get pairwise genetic distance of two populations
    def _get_gd_and_p(self, pops, gd_func):

        pops_haplotype_data = self._get_population_data(pops, self.population_array) # haplotype data of specific populations
        population_data = self.population_array[list(map(lambda x: x in pops, self.population_array[:, self.label_num]))] # population fata of specific populations

        # calculate genetic distance (Fst, Rst or Dst)
        if gd_func == self._rst:
            inds = list(range(len(pops_haplotype_data)))
            distance_matrix = np.zeros(shape=(len(inds), len(inds)))
            for ind1 in inds:
                seq1 = pops_haplotype_data[ind1]
                for ind2 in inds:
                    if ind1 == ind2:
                        distance = 0.0
                    else:
                        seq2 = pops_haplotype_data[ind2]
                        distance = self._pairwise_distance(seq1, seq2, 'sum_of_squared_diff')
                    distance_matrix[ind1, ind2] = distance
            gd = gd_func(pops, population_data, distance_matrix)
        else:
            gd = gd_func(pops, population_data, pops_haplotype_data)

        # permutation to get genetic distance
        permuted_index = list(range(0, len(population_data)))
        permuted_population_data = population_data.copy()
        permuted_gds = []
        for i in range(1000):
            random.shuffle(permuted_index)
            permuted_pops = population_data[permuted_index, self.label_num]
            permuted_population_data[:, self.label_num] = permuted_pops
            if gd_func == self._rst:
                permuted_gd = gd_func(pops, permuted_population_data, distance_matrix)
            else:
                permuted_gd = gd_func(pops, permuted_population_data, pops_haplotype_data)
            permuted_gds.append(permuted_gd)

        # calculate p-value
        p_value = np.sum(np.array(permuted_gds) >= gd) / (1000+1)
        p_value = round(p_value, 5)
        if p_value == 0:
            lowest_p_value = round(1/1001, 5)
            p_value = '<%s' % lowest_p_value

        return gd, p_value

    # get pairwise and global genetic distance of all populations
    def _calculate_gd(self, gd_type):

        # create empty gnetic diastance matrix
        gd_df = pd.DataFrame(index=self.pops, columns=self.pops)
        # function used to calculate genetic distances
        if gd_type == 'fst':
            gd_func = self._fst
        elif gd_type == 'rst':
            gd_func = self._rst
        elif gd_type == 'Ngst':
            gd_func = self._Nei_gst
        elif gd_type == 'NCgst':
            gd_func = self._Nei_Chesser_gst
        elif gd_type == 'Hgst':
            gd_func = self._Hedrick_gst

        # calculate pairwise genetic distances
        pop_num = 0
        for pop1 in self.pops:
            pop_num += 1
            for pop2 in self.pops[pop_num:]:
                pops = [pop1, pop2]
                # calculate pairwise genetic distance and p value
                p_gd, p_gd_p = self._get_gd_and_p(pops, gd_func)
                gd_df.at[pop1, pop2] = p_gd_p
                gd_df.at[pop2, pop1] = p_gd

        gd_df = gd_df.fillna('-')

        # calculate global genetic distances
        g_gd, g_gd_p = self._get_gd_and_p(self.pops, gd_func) # calculate global fst or rst and p value

        return gd_df, g_gd, g_gd_p

    # get df of AMOVA
    def _amova_df(self, amova_pop_data, N, total_level_num, levels):

        # pairwise: df between, df within, df total
        # total: df among, df within, df total
        df = ['']*(total_level_num+2)
        df[:total_level_num] = [len(set(amova_pop_data[:, i])) for i in levels[::-1]]
        df[total_level_num] = N
        df[1: total_level_num+1] = [df[i]-df[i-1] for i in range(1, total_level_num+1)] # df within
        df[0] -= 1 # df between/among
        df[-1] = N - 1 # df total
        df = np.array(df)

        return df

    # get SSD of AMOVA
    def _amova_SSD(self, amova_pop_data, distance_matrix, N, total_level_num, levels):

        SSD = ['']*(total_level_num+2) # SSD between/among, SSD within, SSD total
        SSD[-1] = np.sum(distance_matrix)/float(2*N) # SSD total
        # calculate SSD within
        for level_num, level_name in enumerate(levels[::-1]):
            level_labels = amova_pop_data[:, level_name].tolist()
            labels_num = []
            for i in Counter(level_labels).values():
                labels_num += [i]*i
            level_labels_code = [sorted(list(set(level_labels))).index(i) for i in level_labels]
            SSD[level_num+1] = np.sum((distance_matrix/(2*np.array(labels_num)))[np.equal.outer(level_labels_code, level_labels_code)])
        if total_level_num == 2: # Group in levels
            SSD[1] = SSD[1] - SSD[2]
        SSD[0] = SSD[-1] - sum(SSD[1:-1]) # SSD between/among

        return SSD

    # get coefficient of variance in AMOVA
    def _amova_coef(self, amova_pop_data, N, total_level_num):

        pop_num = len(set(amova_pop_data[:, 1]))
        pop_level_count = np.array(list(Counter(amova_pop_data[:, 1]).values()))
        pop_level_count_square = pop_level_count**2

        if total_level_num == 1:
            # calculate coefficient of variance a
            coef = (N-np.sum(pop_level_count_square)/N) / (pop_num-1)
        elif total_level_num == 2:
            coef = [''] * 3
            group_num = len(set(amova_pop_data[:, 2]))
            # count how many populations belong to a group
            # get each first group-pop of group
            new_array = [tuple(row) for row in amova_pop_data[:, 1:]]
            group_pop = np.unique(new_array, axis=0)[:, 2-1]
            group_count_in_group_pop = list(Counter(group_pop).values()) # the count of group
            # count the sum of population freq in each group
            pop_count_in_group = []
            pop_count_in_group_squre = []
            group_labels = sorted(list(set(amova_pop_data[:, 2])))
            for i in group_labels:
                group_data = amova_pop_data[amova_pop_data[:, 2]==i]
                group_count = np.array(list(Counter(group_data[:, 1]).values()))
                pop_count_in_group.append(np.sum(group_count))
                pop_count_in_group_squre.append(np.sum(group_count**2))
            ratio = np.sum(np.array(pop_count_in_group_squre) / np.array(pop_count_in_group))
            # calculate coefficient of variance
            coef[0] = (N-ratio) / np.sum(np.array(group_count_in_group_pop)-1) # n1
            coef[1] = (ratio-np.sum(pop_level_count_square)/N) / (group_num-1) # n2
            coef[2] = (N-np.sum(np.array(pop_count_in_group)**2/N)) / (group_num-1) # n3
        coef = np.round(np.array(coef), 5)

        return coef

    # get Fst or Rst based on AMOVA analysis
    # get AMOVA table
    def _amova(self, pops, type, global_flag=False):

        if global_flag:
            if self.label == 'Population':
                levels = [self.label_num]
            else:
                levels = list(range(1, np.size(self.population_array, 1)))
            amova_pop_data = self.population_array
        else:
            levels = [self.label_num]
            amova_pop_data = self.population_array[list(map(lambda x: x in pops, self.population_array[:, self.label_num]))]

        amova_inds = amova_pop_data[:, 0]
        amova_haplotype_data = self.haplotype_array[list(map(lambda x: x in amova_inds, self.haplotype_array_with_inds[:, 0]))]
        total_level_num = len(levels) # number of levels: population or group
        N = len(amova_pop_data) # sample size

        # get distance matrix of two populations
        distance_matrix = np.zeros(shape=(len(amova_inds), len(amova_inds)))
        for ind1 in range(len(amova_inds)):
            seq1 = amova_haplotype_data[ind1]
            for ind2 in range(len(amova_inds)):
                if ind1 == ind2:
                    distance = 0
                else:
                    seq2 = amova_haplotype_data[ind2]
                    distance = self._pairwise_distance(seq1, seq2, type)
                distance_matrix[ind1, ind2] = distance

        # calculate df
        df = self._amova_df(amova_pop_data, N, total_level_num, levels)

        # calculate SSD
        SSD = self._amova_SSD(amova_pop_data, distance_matrix, N, total_level_num, levels)

        # calculate MSD
        MSD = SSD / df

        # # calculate coefficient and variance
        coef = self._amova_coef(amova_pop_data, N, total_level_num)
        sigma = ['']*(total_level_num+1)
        if total_level_num == 1:
            sigma[1] = MSD[1] # sigma b
            sigma[0] = (MSD[0]-MSD[1])/coef # sigma a
            statistics = [sigma[0] / (sigma[0]+sigma[1])] # phi_ST
        elif total_level_num == 2:
            sigma[2] = MSD[2] # sigma c
            sigma[1] = (MSD[1]-sigma[2]) / coef[0] # sigma b
            sigma[0] = (MSD[0]-MSD[2]-coef[1]*sigma[1]) / coef[2] # sigma a
            statistics = ['']*3
            statistics[0] = (sigma[0]) / (sigma[0]+sigma[1]+sigma[2]) # phi_CT
            statistics[1] = sigma[1] / (sigma[0]+sigma[1]) # phi_SC
            statistics[2] = (sigma[0]+sigma[1]) / (sigma[0]+sigma[1]+sigma[2]) # phi_ST
        sigma = np.round(np.array(sigma), 5)
        statistics = np.round(np.array(statistics), 5)

        # permutaion to get p value
        p_value = ['']*(total_level_num)
        r_sigma =  np.zeros(shape=(1000, total_level_num)) # sigma number = level number
        if total_level_num == 1:
            for i in range(1000):
                perm_order = np.random.permutation(np.arange(N)) # random permutation order
                r_distance_matrix = distance_matrix[np.ix_(perm_order, perm_order)] # permuted distance matrix
                r_SSD = self._amova_SSD(amova_pop_data, r_distance_matrix, N, total_level_num, levels)
                r_MSD = r_SSD / df
                r_sigma[i, 0] = (r_MSD[0]-r_MSD[1])/coef
        elif total_level_num == 2:
            group_label_count = list(Counter(amova_pop_data[:, 2]).values())
            end = 0
            label_num = []
            for i, j in enumerate(group_label_count):
                end += j
                start = sum(group_label_count[:i])
                label_num.append([start, end])
            for i in range(1000):
                perm_order = np.array([])
                for num in label_num:
                    perm_order = np.append(perm_order, np.random.permutation(np.arange(num[0], num[1])))
                perm_order = perm_order.astype('int64')
                r_distance_matrix = distance_matrix[np.ix_(perm_order, perm_order)]
                r_SSD = self._amova_SSD(amova_pop_data, r_distance_matrix, N, total_level_num, levels)
                r_MSD = r_SSD / df
                r_sigma[i, 1] = (r_MSD[1]-r_MSD[2]) / coef[0]
            pop_label_count = list(Counter(amova_pop_data[:, 1]).values())
            new_array = [tuple(row) for row in amova_pop_data[:, levels]]
            pop_group = np.unique(new_array, axis=0)[:, 2-1]
            r_amova_pop_data = amova_pop_data.copy()
            for i in range(1000):
                random.shuffle(pop_group)
                r_amova_pop_data[:, 2] = list(chain.from_iterable(repeat(group, count) for group, count in zip(pop_group, pop_label_count)))
                r_SSD = self._amova_SSD(r_amova_pop_data, r_distance_matrix, N, total_level_num, levels)
                r_df = self._amova_df(r_amova_pop_data, N, total_level_num, levels)
                r_MSD = r_SSD / r_df
                r_coef = self._amova_coef(r_amova_pop_data, N, total_level_num)
                r_sigma[i, 0] = (r_MSD[0]-r_MSD[2]-r_coef[1]*((r_MSD[1]-r_MSD[2])/r_coef[0])) / r_coef[2]
        for i in range(total_level_num):
            p_value[i] = np.round(np.sum(r_sigma[:, i] >= sigma[i]) / (1000+1), 5)

        if 0 in p_value:
            lowest_p_value = round(1/1001, 5)
            p_value = ['<%s' % lowest_p_value if i == 0 else i for i in p_value]

        if global_flag:
            return SSD, df, MSD, sigma, p_value, statistics, coef
        else:
            return statistics[0], p_value[0]

    # calculate pairwise Fst or Rst of populations
    def _calculate_amova(self, stat_type):

        if stat_type == 'fst':
            type = 'num_of_diff_alleles'
        elif stat_type == 'rst':
            type = 'sum_of_squared_diff'

        stat_df = pd.DataFrame(index=self.pops, columns=self.pops)
        pop_num = 0
        for pop1 in self.pops:
            pop_num += 1
            for pop2 in self.pops[pop_num:]:
                pops = [pop1, pop2]
                stat, p = self._amova(pops, type)
                stat_df.at[pop1, pop2] = p
                stat_df.at[pop2, pop1] = stat

        stat_df = stat_df.fillna('-')

        # Global AMOVA
        amova_items = self._amova(self.pops, type, True) # SSD, Df, MSD, sigma, statistics, p_value
        # AMOVA table
        if self.label == 'Population':
            amova_df_index = ['Among_populations', 'Within_populations', 'Total']
            variance_df_index = ['a', 'b']
            phi_df_columns = ['Phi_ST']
            coef_df_columns = ['n1']
        elif self.label == 'Group':
            amova_df_index = ['Among_groups', 'Among_populations/Within_Groups', 'Within_populations', 'Total']
            variance_df_index = ['a', 'b', 'c']
            phi_df_columns = ['Phi_CT', 'Phi_SC', 'Phi_ST']
            coef_df_columns = ['n1', 'n2', 'n3']
        amova_df = pd.DataFrame(index=amova_df_index, columns=['SSD', 'Df', 'MSD'])
        variance_df = pd.DataFrame(index=variance_df_index, columns=['Variance', 'p_value'])
        phi_df = pd.DataFrame(index=['Phi_statistic'], columns=phi_df_columns)
        coef_df = pd.DataFrame(index=['Variance_coefficients'], columns=coef_df_columns)
        for i, j in zip(amova_items[:-2], ['SSD', 'Df', 'MSD', 'Variance', 'p_value']):
            if j in amova_df.columns:
                tmp_df = amova_df
            else:
                tmp_df = variance_df
            if len(i) != tmp_df.index.size:
                diff = tmp_df.index.size - len(i)
                tmp_df[j] = np.append(i, np.array(['']*diff))
            else:
                tmp_df[j] = i

        phi_df.loc['Phi_statistic'] = amova_items[-2]
        coef_df.loc['Variance_coefficients'] = amova_items[-1]

        return stat_df, amova_df, variance_df, phi_df, coef_df

    # output calaulation result
    def calculate_and_output(self, arguments, path):

        if arguments.hd:
            hd_df = self._calculate_haplotype_diversity()
            hd_df.to_csv(path+'.HD.%s.txt' % self.label, sep='\t')
            self.logger.info('[Y-LineageTracker] Calculate haplotype diversity finished')
        if arguments.mpd:
            mpd_df = self._calculate_mean_pairwise_differences()
            mpd_df.to_csv(path+'.MPD.%s.txt' % self.label, sep='\t')
            self.logger.info('[Y-LineageTracker] Calculate mean pairwise differences on %s level finished' % self.label)
        if arguments.gd:
            self.logger.info('[Y-LineageTracker] Calculating pairwise genetic distance on %s level and p-value from permutation……' % self.label)
            gd_df, global_gd, global_gd_p = self._calculate_gd(arguments.gd)
            gd_out = open(path+'.GD.%s.%s.txt' % (arguments.gd, self.label), 'w')
            gd_out.write(gd_df.to_string()+'\n\n')
            gd_out.write('Global %s: %f\n' % (arguments.gd, global_gd))
            gd_out.write('p-value of Global %s: %s' % (arguments.gd, global_gd_p))
            gd_out.close()
            self.logger.info('[Y-LineageTracker] Calculation of pairwise genetic distance and p-value on %s level finished' % self.label)
        if arguments.amova:
            self.logger.info('[Y-LineageTracker] Calculating pairwise %s based on AMOVA of %s level……' % (arguments.amova, self.label))
            stat_df, amova_df, variance_df, phi_df, coef_df = self._calculate_amova(arguments.amova)
            stat_df.to_csv(path+'.AMOVA.%s.%s.txt' % (arguments.amova, self.label), sep='\t')
            self.logger.info('[Y-LineageTracker] Calculating pairwise %s based on AMOVA of %s level finished' % (arguments.amova, self.label))
            amova_df_out = open(path+'.AMOVA.%s.txt' % (self.label), 'w')
            amova_df_out.write(amova_df.to_string()+'\n\n'+variance_df.to_string()+'\n\n'+phi_df.to_string()+'\n\n'+coef_df.to_string())
            amova_df_out.close()
            self.logger.info('[Y-LineageTracker] AMOVA on %s level finished' % self.label)


# predict Y-haplogroup from Y-STRs
def predict_haplogroup(haplotype_data, path):

    logger = logging.getLogger()
    skip_letters = config.get('Statistics', 'SkipLetters').split(',')
    header = haplotype_data.columns.tolist()
    haplotype_array = np.array(haplotype_data)

    # get Y-STR prediction panel
    from FilesIO import CommonData
    common_data = CommonData()
    predict_panel = common_data.read_str_prediction_panel()
    array_columns = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'L', 'M', 'N', 'O', 'Q', 'R', 'S', 'T']

    to_be_imputed = 'Haplogroup'
    predict_panel[to_be_imputed] = predict_panel[to_be_imputed].map(lambda x: x[0])

    code, code2 = 0, 0
    str_list, str_list2 = [], []
    sample_list = []
    array_list = []
    str_panel = {}
    used = {}
    str_info = {}
    # get samples and their Y-STRs
    for n, sample in enumerate(haplotype_data.index):
        STR = {}
        for m, STR_name in enumerate(header):
            STR_num = haplotype_array[n, m]
            if STR_num not in skip_letters:
                try:
                    STR[STR_name] = int(STR_num)
                except ValueError:
                    logger.wrning('[Y-LineageTracker] [Warning] The STR input of %s does not meet the requirements' % sample)

        # prediction panel for each ind
        STRs = list(STR.keys())
        if STRs in str_list:
            code = str_list.index(STRs)
            used_STRs = used[code]
            ind_predict_panel = str_panel[code]
        else:
            na_count = predict_panel[STRs].isnull().sum(axis=0)
            used_STRs = na_count[na_count<1000].index.tolist()
            if len(used_STRs) < 6:
                ind_predict_panel = None
            else:
                ind_predict_panel = predict_panel[[to_be_imputed]+used_STRs]
            used[code] = used_STRs
            str_panel[code] = ind_predict_panel
            str_list.append(STRs)
            code += 1
        if len(used_STRs) < 6:
            logger.warning('[Y-LineageTracker] [Warning] The number of Y-STR loci is too small for sample %s' % sample)
            continue

        freq_list = []

        # calculate frequency
        ind_predict_hap = np.array(ind_predict_panel[to_be_imputed])
        STR_loci = ind_predict_panel.drop(to_be_imputed, axis=1).columns.tolist()
        ind_predict_panel = np.array(ind_predict_panel.drop(to_be_imputed, axis=1))

        ind_ref_panel = np.empty(shape=(0, np.size(ind_predict_panel, 1)))
        ind_ref_hap = []
        for i in sorted(list(set(ind_predict_hap))):
            sub = ind_predict_panel[ind_predict_hap==i]
            if len(sub) > 500:
                idx = np.random.randint(0, len(sub), 500)
                sub = sub[idx]
            ind_ref_panel = np.concatenate([ind_ref_panel, sub])
            ind_ref_hap.extend([i]*len(sub))
        ind_ref_hap = np.array(ind_ref_hap)
        to_be_imputed_freq = np.array(list(zip(*sorted(Counter(ind_ref_hap).items())))[1]).astype(float)

        # sample size control
        panel = {}
        array_location = {}
        for i, locus in enumerate(STR_loci):
            array_index = sorted(np.unique(ind_ref_panel[:, i]))
            array = np.zeros(shape=(len(array_index), len(array_columns)))
            for j, hap in enumerate(array_columns):
                col_data = ind_ref_panel[list(map(lambda x: x == hap, ind_ref_hap))][:, i]
                col_freq = Counter(col_data)
                col = [col_freq[i]/len(col_data) if i in col_freq.keys() else 0 for i in array_index]
                array[:, j] = col
            panel[locus] = array
            array_location[locus] = array_index

        # calculate the probability of each main haplogroup
        Pr_up = to_be_imputed_freq
        for i in STR.keys():
            if i in used_STRs:
                table_index = array_location[i].index(STR[i])
                Pr_up *= panel[i][table_index]
        D = sum(Pr_up)
        Pr = np.round(Pr_up/D, 4)
        array_list.append(Pr)
        sample_list.append(sample)

        percent = float(n + 1) * 100 / float(len(haplotype_data))
        sys.stdout.write('[Y-LineageTracker] Predicting haplogroups from Y-STR data…… %.2f' % percent)
        sys.stdout.write('%\r')
        sys.stdout.flush()

    # output result
    predict_df = pd.DataFrame(index=sample_list,
                              columns=sorted(list(set(predict_panel['Haplogroup'].map(lambda x: x[0])))),
                              data=np.array(array_list))
    predict_df.to_csv(path+'.Prediction.txt', sep='\t')

    print()
    logger.info('[Y-LineageTracker] Haplogroup prediction finished')


def main():

    start = time.perf_counter()
    arguments = stat_parser()

    from FilesIO import check_seq_matrix_input, check_hap_input, set_haplotype_index, check_population, check_overlap, check_pop_header, time_count
    # check data input
    input, format = check_seq_matrix_input(arguments.seq, arguments.format, arguments.matrix)
    haplotype_data = check_hap_input(input, 'haplotype', format)

    # check commands
    arguments.amova = check_amova(arguments.amova, sys.argv)
    check_commands(arguments)

    from FilesIO import get_out_path
    # get path of output files
    path = get_out_path(arguments.output, 'Stat')
    log_file = path + '.StatLog.log'

    set_log(arguments, log_file)

    # check population data
    population_data = check_population(arguments.population)
    population_data = check_pop_header(population_data)
    haplotype_data, population_data = check_overlap(haplotype_data, population_data)
    haplotype_data = set_haplotype_index(haplotype_data, format)

    # calculate statistics
    label = 'Population'
    stat = StatisticsAnalysis(label, haplotype_data, population_data)
    stat.calculate_and_output(arguments, path)
    if 'Group' in population_data.columns:
        label = 'Group'
        stat = StatisticsAnalysis(label, haplotype_data, population_data)
        stat.calculate_and_output(arguments, path)

    # predict haplogroup from Y-STR data
    if arguments.predict and arguments.matrix:
        predict_haplogroup(haplotype_data, path)

    time_count(start)


if __name__ == '__main__':
    main()
