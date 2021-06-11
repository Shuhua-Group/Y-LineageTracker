import os
import sys
import time
import logging
import argparse
import matplotlib; matplotlib.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from itertools import filterfalse
from collections import Counter
from matplotlib import font_manager
from GetConfig import getConfig

try:
    from adjustText import adjust_text
except ImportError:
    print('Module [adjustText] is required, please install it by command line: pip install adjustText')
    sys.exit()


'''
Description:
This module is used for clustering analysis.
Clustering method can be applied by --method option.
The program will take haplogroup file and population file as input files.
Then cluster populations by PCA or MDS method and visualize results by the plot.
'''


config = getConfig() # get config


def cluster_parser():

    parser = argparse.ArgumentParser('cluster', description='(c) Y-LineageTracker: NRY haplogroup clustering analysis')
    # function used for clustering analysis
    parser.add_argument('cluster',
                        help='Pefrorm clustering analysis for NRY haplogroups.')
    # required, haplogroup file
    parser.add_argument('--hg',
                        required=True,
                        type=str,
                        action='store',
                        help='hg: A file containing sample ID and haplogroup of each individual.')
    # required, population file
    parser.add_argument('-p', '--population',
                        required=True,
                        type=str,
                        action='store',
                        help='population: A file containing sample ID and population information of each individual.')
    # optional, clustering algorithm
    parser.add_argument('--method',
                        required=False,
                        type=str,
                        action='store',
                        default='pca',
                        choices=['pca', 'mds'],
                        help='method: The algorithm used for clustering analysis (PCA or MDS).')
    # optional, simpification mode for haplogroups in same main trunk
    parser.add_argument('--level',
                        required=False,
                        type=str,
                        const='auto',
                        action='store',
                        nargs='?',
                        help='Level: the haplogroup resolution level of NRY haplogroups for clustering,\
                              which means the degree of similarity among haplogroups in the same main trunk.\
                              The parameter [auto] will set an optimal level to simplify NRY haplogroups resolution.')
    # optional, output haplogroup frequence of populations
    parser.add_argument('--freq',
                        required=False,
                        action='store_true',
                        help='frequence: Output the haplogroup frequency of each population.')
    # optional, the prefix of output
    parser.add_argument('-o', '--output',
                        required=False,
                        type=str,
                        action='store',
                        help='output: The prefix of output files.')


    args = parser.parse_args()

    return args


# print program information and write to log file
def set_log(log_file, level, args_log):

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

    log_info = ['[Y-LineageTracker] [Cluster]',
                '[Y-LineageTracker] Run Date: ' + time.asctime(time.localtime(time.time())),
                '[Y-LineageTracker] Input File: %s' % args_log.hg,
                '[Y-LineageTracker] Clustering Method: %s' % args_log.method,
                '[Y-LineageTracker] Haplogroup Level: %s' % level]

    if args_log.freq:
        log_info.append('[Y-LineageTracker] Output haplogroup frequence')

    print('\n')
    for i in log_info:
        logger.info(i)


# names of output files
def check_cluster_output(output, method, freq):

    from FilesIO import get_out_path
    path = get_out_path(output, 'Cluster')
    level_file = path + '.level.hg'
    if method == 'pca':
        fig_file = path + '.PCA.fig.pdf'
        vect_file = path + '.PCA.eigenvectors.txt'
        val_file = path + '.PCA.eigenval.txt'
        output_set = [fig_file, vect_file , val_file, level_file]
        log_file = path + '.ClusterLog.PCA.log'
    elif method == 'mds':
        fig_file = path + '.MDS.fig.pdf'
        fst_file = path + '.MDS.fst.txt'
        embedding_file = path + '.MDS.embedding.txt'
        output_set = [fig_file, fst_file, embedding_file, level_file]
        log_file = path + '.ClusterLog.MDS.log'
    if freq:
        freq_file = path + '.freq.txt'
        output_set.append(freq_file)

    return output_set, log_file


# function of haplogroup simplification
def hap_simplify(haplogroup, level):

    if haplogroup.startswith('A'):
        if haplogroup.startswith('A1'):
            return haplogroup[:level+1]
        elif haplogroup.startswith('A0'):
            zero_num = haplogroup.count('0')
            if 'T' in haplogroup:
                return haplogroup[:level+zero_num+1]
            else:
                return haplogroup[:level+zero_num]
    else:
        if (haplogroup.isalpha() and haplogroup.isupper()):
            return haplogroup
        elif haplogroup[:-1].isupper and haplogroup[:-1].isalpha() and len(haplogroup) > 2:
            return haplogroup
        else:
            return haplogroup[:level]

# check level parameter
def check_level(level):

    if level:
        if level in ['auto', 'original', 'min']:
            return level
        else:
            try:
                level = int(level)
                if level >= 1:
                    return level
                else:
                    raise argparse.ArgumentTypeError('parameter [level] should be an integer or auto')
            except ValueError:
                raise argparse.ArgumentTypeError('parameter [level] should be an integer or auto')
    else:
        level = 'original'
        return level


class ClusterAnalysis(object):
    '''
    This class is used to perform main clustering analysis
    Four main steps are included in this class:
    1. read haplogroup file and simplify NRY haplogroups
    2. read population file and merge to haplogroup data
    3. perform clustering analysis
    4. output clustering result with figure
    '''

    def __init__(self):

        self.logger = logging.getLogger()

    # calculate Fst
    # Weir and Hill. 2002. Estimating F-statistics. Ann.Rev.Gen.
    def _calculate_pairwise_fst(self, pop1_count, pop2_count):

        sample_size1 = (np.sum(pop1_count)/pop1_count.size)*2
        sample_size2 = (np.sum(pop2_count)/pop2_count.size)*2

        pop1_frequency = np.array(pop1_count) / sample_size1
        pop2_frequency = np.array(pop2_count) / sample_size2

        mean_freq = (np.array(pop1_count)+np.array(pop2_count)) / (sample_size1+sample_size2)

        MSP1 = sample_size1*((pop1_frequency-mean_freq)**2)
        MSP2 = sample_size2*((pop2_frequency-mean_freq)**2)
        MSP = MSP1 + MSP2

        MSG1 = float(1/((sample_size1-1)+(sample_size2-1)))
        MSG2 = (np.array(pop1_count) * (1-pop1_frequency)) + (np.array(pop2_count) * (1-pop2_frequency))
        MSG = MSG1 * MSG2

        nc = sample_size1+sample_size2 - (sample_size1**2+sample_size2**2)/(sample_size1+sample_size2)

        theta1 = np.sum(MSP-MSG)
        theta2 = np.sum(MSP+(nc-1)*MSG)
        if theta2 == 0:
            theta = 0
        else:
            theta = theta1 / theta2

        return theta

    # calculate pairwise Fst
    def _calculate_MDS_fst(self, population_count, fst_output):

        pops = population_count.index.tolist()
        fst_df = pd.DataFrame(index=pops, columns=pops)
        pop_num = 0
        for pop1 in pops:
            for pop2 in pops[pop_num:]:
                if pop1 == pop2:
                    fst_df.at[pop2, pop1] = 0
                else:
                    pop1_count = np.empty(shape=(2, population_count.columns.size))
                    pop1_count[0] = np.array(population_count.loc[pop1].tolist()).reshape(1, -1)
                    pop1_count[1] = (np.sum(pop1_count[0])-np.array(population_count.loc[pop1].tolist())).reshape(1, -1)
                    pop2_count = np.empty(shape=(2, population_count.columns.size))
                    pop2_count[0] = np.array(population_count.loc[pop2].tolist()).reshape(1, -1)
                    pop2_count[1] = (np.sum(pop2_count[0])-np.array(population_count.loc[pop2].tolist())).reshape(1, -1)
                    fst = self._calculate_pairwise_fst(pop1_count, pop2_count)
                    fst_df.at[pop2, pop1] = fst
                    fst_df.at[pop1, pop2] = fst
            pop_num += 1

        fst_df.to_csv(fst_output, sep='\t', index=True)
        fst_array = np.array(fst_df)

        return fst_array

    # simplify haplogroups according to level
    def classify_hap_data(self, level, hap_data, level_file):

        if level != 'original':
            # simplify to most common haplogroup
            if level == 'auto':
                initial = set(hap_data['Haplogroup'].map(lambda x: x[0]))
                common_list = []
                for i in initial:
                    initial_group = sorted(list(filter(lambda x: x[0] == i, hap_data['Haplogroup'])), key=lambda x: len(x))
                    while True:
                        most_common = Counter(initial_group).most_common(1)[0][0]
                        initial_group = [most_common if x.startswith(most_common) and len(x) > len(most_common) else x for x in initial_group]
                        hap_num = len(most_common)-1
                        for j in range(len(most_common)-1):
                            hap = most_common[0:hap_num]
                            if hap in initial_group:
                                initial_group = [hap if x == most_common else x for x in initial_group]
                                most_common = hap
                            hap_num -= 1
                        common_list.append(most_common)
                        initial_group = list(filterfalse(lambda x: x.startswith(most_common), initial_group))
                        if len(initial_group) == 0:
                            break
                common_list = sorted(common_list, key=lambda x: len(x), reverse=True)
                for i in range(len(hap_data['Haplogroup'])):
                    for j in common_list:
                        if hap_data.at[i, 'Haplogroup'].startswith(j):
                            hap_data.at[i, 'Haplogroup'] = j
            # simplify to specific level
            else:
                if level == 'min':
                    main_trunk = config.get('HaplogroupTree', 'CommonTrunk').split(',') + config.get('HaplogroupTree', 'UpperTrunk').split(',')
                    level = min(map(lambda x: len(x), fiter(lambda x: x not in main_trunk, hap_data['Haplogroup'])))
                hap_data['Haplogroup'] = hap_data['Haplogroup'].apply(hap_simplify, level=level)
            hap_data.to_csv(level_file, sep='\t', index=False)


        return hap_data

    # read population data and plot clustering figure
    def output_cluster(self, method, freq, output_set, hap_data, population_data):

        #############################READ POPULATION DATA#############################

        # merge population data for frequence construction
        hap_data = pd.merge(hap_data, population_data, on='SampleID')

        # count number of populations and sort them in order
        num_list = []
        num = 0
        for i in sorted(list(set(population_data['Population'])), key=population_data.drop('SampleID', axis=1).drop_duplicates()['Population'].tolist().index):
            num_list.append((num, num+population_data['Population'].tolist().count(i)))
            num = num + population_data['Population'].tolist().count(i)
        population_data = population_data.drop('SampleID', axis=1).drop_duplicates()

        # sort according to Group for subsequent plot
        if 'Group' in population_data.columns:
            population_data = population_data.sort_values(by='Group')

        # haplogorup frequence of populations
        population_freq = pd.DataFrame(index=sorted(list(set(hap_data['Population'])),
                                       key=population_data['Population'].tolist().index),
                                       columns=sorted(list(set(hap_data['Haplogroup']))),
                                       dtype=np.float)
        hap_series = pd.Series(hap_data['Haplogroup'].tolist(), index=[hap_data['Population']])
        if method == 'mds':
            population_count = pd.DataFrame(index=sorted(list(set(hap_data['Population'])),
                                            key=population_data['Population'].tolist().index),
                                            columns=sorted(list(set(hap_data['Haplogroup']))),
                                            dtype=np.int)
        # count haplogroup frequence
        for i in population_freq.index:
            hap_list = hap_series[i].value_counts().index.tolist()
            freq_list = hap_series[i].value_counts().tolist()
            sum_num = sum(freq_list)
            for j in range(len(hap_list)):
                population_freq.at[i, hap_list[j]] = freq_list[j] / sum_num
                if method == 'mds':
                    population_count.at[i, hap_list[j]] = freq_list[j]
        population_freq = population_freq .fillna(0)
        if method == 'mds':
            population_count = population_count.fillna(0)

        #############################PERFORM CLUSTERING ALGORITHM#############################

        matrix = population_freq.to_numpy() # convert haplogroup frequency dataframe to martix
        if method == 'pca': # PCA analysis
            from sklearn.decomposition import PCA
            pca = PCA(copy=True, whiten=False)
            coordinates = pca.fit_transform(matrix)
            # output PCA results
            output_vect_file = output_set[1]
            output_vect_data = pd.DataFrame(coordinates, index=sorted(list(set(hap_data['Population']))))
            output_vect_data.to_csv(output_vect_file , sep='\t', index=True)
            output_val_file = output_set[2]
            output_val_data = pd.DataFrame(pca.explained_variance_)
            output_val_data.to_csv(output_val_file, index=False, header=None)
        elif method == 'mds': # MDS analysis
            from sklearn.manifold import MDS
            mds_matrix = self._calculate_MDS_fst(population_count, output_set[1])
            mds = MDS(dissimilarity='precomputed')
            coordinates = mds.fit_transform(mds_matrix)
            # output MDS result
            output_embedding_file = output_set[2]
            output_embedding_data = pd.DataFrame(coordinates, index=sorted(list(set(hap_data['Population']))))
            output_embedding_data.to_csv(output_embedding_file , sep='\t', index=True)
        # output frequence of haplogroup in populations
        if freq:
            output_freq_file = output_set[-1]
            population_freq.sort_index().to_csv(output_freq_file, sep='\t', index=True)

        #############################START PLOT#############################

        # defualt color and shape in cluster if no Group columns in population file
        color_list = ['slategrey']*population_data.index.size
        shape_list = ['o']*population_data.index.size

        # set color according to Group labels
        if 'Group' in population_data.columns:
            label_list = population_data['Group'].tolist()
            # if with color columns
            if 'Color' in population_data.columns:
                color_list = population_data['Color'].tolist()
            # get cmap colors in matplotlib
            else:
                from FilesIO import set_color_num
                color_num = len(set(label_list))
                colors = set_color_num(color_num)
                label_count = 0
                color_list = []
                colored_label = []
                label_dict = {}
                for i in label_list:
                    if i not in colored_label:
                        label_dict[i] = colors[label_count]
                        colored_label.append(i)
                        color_list.append(label_dict[i])
                        label_count += 1
                    else:
                        color_list.append(label_dict[i])
        # No Group columns
        elif 'Color' in population_data.columns and 'Group' not in population_data.columns:
            color_list = population_data['Color'].tolist()
        if 'Shape' in population_data.columns:
            shape_list = population_data['Shape'].tolist()

        self.logger.info('[Y-LineageTracker] Start clustering analysis')

        # parameters for plot
        import warnings
        warnings.filterwarnings('ignore')
        font_files = font_manager.findSystemFonts(fontpaths=[os.path.split(os.path.realpath(__file__))[0] + '/sans-serif'])
        font_list = font_manager.createFontList(font_files)
        font_manager.fontManager.ttflist.extend(font_list)
        plt.rcParams['font.family']= 'sans-serif'
        plt.rcParams['font.sans-serif'] = ['Arial'] # font
        plt.rcParams['axes.unicode_minus'] = False
        figure, ax = plt.subplots(figsize = (12, 12)) # figure size

        # x and y lines of zero
        ax.axvline(x=0, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)
        ax.axhline(y=0, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)

        # x and y labels
        if method == 'pca':
            plt.xlabel('PC1  %.2f%%'% (pca.explained_variance_ratio_[0]*100), fontsize=25)
            plt.ylabel('PC2  %.2f%%'% (pca.explained_variance_ratio_[1]*100), fontsize=25)
        elif method == 'mds':
            plt.xlabel('Dimension1', fontsize=25)
            plt.ylabel('Dimension2', fontsize=25)

        plt.xticks(size=20)
        plt.yticks(size=20)

        # traverse populations to plot scatterplot
        if 'label_list' in vars():
            plot_list = []
            for i in range(population_data.index.size):
                plot = list(zip(color_list, shape_list, label_list))[i]
                if plot in plot_list:
                    plt.plot(coordinates[i, 0], coordinates[i, 1], shape_list[i],
                             color=color_list[i],
                             mec='black',
                             markersize=15,
                             alpha=0.9)
                else:
                    plt.plot(coordinates[i, 0], coordinates[i, 1], shape_list[i],
                             color=color_list[i],
                             mec='black',
                             label=label_list[i],
                             markersize=15,
                             alpha=0.9)
                plot_list.append(plot)
            plt.legend(frameon=True, fontsize=15, markerscale=1, shadow=True, edgecolor='dimgrey', loc='best')
        else:
            for i in range(population_data.index.size):
                plt.plot(coordinates[num_list[i][0]:num_list[i][1], 0], coordinates[num_list[i][0]:num_list[i][1], 1], shape_list[i],
                         color=color_list[i],
                         mec='black',
                         markersize=15,
                         alpha=0.7)

        # plot population with labels
        if population_freq.index.size < 40:
            texts = [plt.text(coordinates[i, 0], coordinates[i, 1], s=population_freq.index.tolist()[i], size=15) for i in range(len(matrix))]
            import warnings
            warnings.filterwarnings('ignore')
            adjust_text(texts)
        else:
            self.logger.warning('[Y-LineageTracker] Too much population in the population file')

        # output figure
        output_fig_file = output_set[0]
        plt.savefig(output_fig_file, dpi=100, bbox_inches='tight')
        plt.close()

        self.logger.info('[Y-LineageTracker] Clustering analysis finished')


def main():

    start = time.perf_counter()
    arguments = cluster_parser()

    # check haplogroup level
    global level
    level = check_level(arguments.level)

    # set of output files
    output_set, log_file = check_cluster_output(arguments.output, arguments.method, arguments.freq)

    # set log file
    set_log(log_file, level, arguments)

    from FilesIO import check_hap_input, check_population, check_overlap, time_count
    # check haplogroup data
    hap_data = check_hap_input(arguments.hg, 'haplogroup')

    # check population data
    population_data = check_population(arguments.population)

    # check overlap of haplogroup data and population data
    hap_data, population_data = check_overlap(hap_data, population_data)

    # start cluserting analysis
    cluster = ClusterAnalysis()
    hap_data = cluster.classify_hap_data(level, hap_data, output_set[3])
    cluster.output_cluster(arguments.method, arguments.freq, output_set, hap_data, population_data)

    time_count(start)


if __name__ == '__main__':
    main()
