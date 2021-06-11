import os
os.environ['QT_QPA_PLATFORM']='offscreen'
import gc
import sys
import time
import logging
import argparse
import matplotlib.colors
import numpy as np
import pandas as pd

from Bio import Phylo
from itertools import filterfalse
from ete3 import Tree
from GetConfig import getConfig


'''
Description:
This module is used for phylogeny analysis.
It is recommended to input both haplogroup file and sequence alignment file to construct a bifurcating tree.
If only the haplogroup file is given, the phylogenetic tree may contain polytomies (a node with more than two children).
For polytomies, if sequence alignment is provided, the program will apply alignment method to construct a sub-tree for nodes.
'''


config = getConfig()


def phylo_parser():

    parser = argparse.ArgumentParser('phylo', description='(c) Y-LineageTracker: Phylogeny analysis')
    # function used for phylogeny analysis
    parser.add_argument('phylo',
                        help='Phylogenetic analysis for NRY haplogroups.')
    # required, haplogroup file
    parser.add_argument('--hg',
                        required=True,
                        type=str,
                        action='store',
                        help='hg: A file containing sample ID and haplogroup of each individual.')
    # optional, sequence alignment file
    parser.add_argument('--seq',
                        required=False,
                        type=str,
                        action='store',
                        help='seq: Sequence alignment file of input samples.')
    # optional, format of sequence alignment
    parser.add_argument('--seq-format',
                        required=False,
                        type=str,
                        dest='format',
                        choices=['fasta', 'phylip', 'nexus', 'meg', 'vcf'],
                        help='seq-format: The file format of sequence alignment file. The default is fasta. ')
    # optional, population file
    parser.add_argument('-p', '--population',
                        required=False,
                        type=str,
                        action='store',
                        help='population: A file containing sample ID and population information of each individual')
    # optional, method of calculating distance from sequence alignment
    parser.add_argument('--align-method',
                        required=False,
                        type=str,
                        dest='align',
                        action='store',
                        default='mp',
                        choices=['ibs', 'upgma', 'mp'],
                        help='align-method: The method of constructing the bifurcating tree from multiple sequences. The default is mp.')
    # optional, the output format in text of tree
    parser.add_argument('--tree-format',
                        required=False,
                        type=str,
                        dest='tree',
                        action='store',
                        default='newick',
                        choices=['newick', 'nexus', 'phyloxml', 'nexml'],
                        help='--tree-format: The file format of output phylogenetic tree. The default is newick. ')
    # optional, the layout of tree in figure
    parser.add_argument('--layout',
                        required=False,
                        type=str,
                        action='store',
                        default='rectangular',
                        choices=['rectangular', 'circular'],
                        help='layout: The layout of the phylogenetic tree figure. The default is rectangular.')
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

    log_info = ['[Y-LineageTracker] [Phylo]',
                '[Y-LineageTracker] Run Date: ' + time.asctime(time.localtime(time.time())),
                '[Y-LineageTracker] Haplogroup File: %s' % args_log.hg,
                '[Y-LineageTracker] Seq File: %s' % args_log.seq,
                '[Y-LineageTracker] Seq Format: %s' % args_log.format,
                '[Y-LineageTracker] Population File: %s' % args_log.population,
                '[Y-LineageTracker] Tree Layout: %s' % args_log.layout,
                '[Y-LineageTracker] Output Tree Format: %s' % args_log.tree]

    print('\n')
    for i in log_info:
        logger.info(i)


# names of output files
def check_phylo_output(tree_format, output):

    from FilesIO import get_out_path
    path = get_out_path(output, 'Phylo')
    fig_file = path + '.phylo'
    if 'xml' in tree_format:
        tree_file = path + '.tree.xml'
    else:
        if tree_format == 'newick':
            tree_file = path + '.tree.nwk'
        elif tree_format == 'nexus':
            tree_file = path + '.tree.nex'
    log_file = path + '.PhyloLog.log'
    output_set = [fig_file, tree_file, log_file]

    return output_set


class PhyloHaplogroup(object):
    '''
    This class is used to perform main phylogeny analysis
    It can:
    1. construct a haplogroup tree form haplogroup file
    2. construct sub-trees for internal polytomies
    3. output tree with a figure
    Population file is recommended for tree plot to construct a bifurcating tree
    Output tree file in newick format can be used as one of input file in time estimation
    '''

    def __init__(self, hap_data, output_set):

        self.logger = logging.getLogger()
        hap_data['Haplogroup'] = hap_data['Haplogroup'].map(lambda x: x.strip('~'))
        self.hap_data = hap_data
        self.output_set = output_set
        self.common_trunk = config.get('HaplogroupTree', 'CommonTrunk').split(',')

    # read sequence alignment file
    def read_seq_file(self, seq_file, format):

        from ProcessData import ConvertData
        from Bio import AlignIO

        # get mutiple alignment data
        if format == 'vcf':
            convert = ConvertData.ConvertData(seq_file)
            msa = convert.data_convert('vcf', 'msa')
        elif format == 'phylip':
            msa = AlignIO.read(open(seq_file), 'phylip-relaxed')
        else:
            msa = AlignIO.read(open(seq_file), format)

        # get ids of msa
        ids = [i.id for i in msa]
        seqs = [i.seq for i in msa]

        return ids, seqs

    # get base tree from haplogroup file
    def get_base_tree(self, used_samples=None):

        self.logger.info('[Y-LineageTracker] Constructing phylogenetic tree...')

        # preprocess data
        hap_df = self.hap_data.sort_values(by='Haplogroup')
        if used_samples:
            hap_df = hap_df[hap_df['SampleID'].map(lambda x: x in used_samples)]
        hap_phylo = hap_df.drop_duplicates(subset='Haplogroup', keep='first', inplace=False)
        hap_inds = hap_df.set_index(['Haplogroup'], drop=True)

        # prepare variables
        new_tree_info = []
        filter_trunk = []
        from FilesIO import CommonData
        common_data = CommonData()
        tree_info = common_data.read_tree_info()

        # traverse tree to build a pylogenetic tree based on haplogroup
        for i in tree_info:
            hap_name = i.split('\t')[-1]
            if hap_name in self.common_trunk:
                new_tree_info.append(i)
                continue
            level = i.count('\t')
            if len(list(filterfalse(lambda x: not x.startswith(hap_name), hap_phylo['Haplogroup']))) == 0:
                new_tree_info.append(i)
                filter_trunk.append(hap_name)
                continue

            hap_point_list = list(filterfalse(lambda x: x == hap_name or not x.startswith(hap_name), hap_phylo['Haplogroup']))
            hap_point_list = [i for i in hap_point_list if i not in self.common_trunk]
            hap_dict = {}
            hap_class = []
            if len(hap_point_list) != 0 and len(set(hap_point_list)&set([i.split('\t')[-1] for i in tree_info])) == 0:
                num = 0
                hap_count = len((hap_point_list))+1
                iter_time = len(''.join(hap_point_list))
                if len(hap_point_list) == 1:
                    new_tree_info.append(i)
                    new_tree_info.append('-\t' * (level+1) + hap_point_list[0])
                    continue
                for j in range(iter_time):
                    if num+1 < len(hap_point_list):
                        length = len(hap_name)
                        if len(hap_point_list[num]) == length:
                            num += 1
                            hap_name = i.split('\t')[-1]
                            hap_count = len(hap_dict[hap_name])
                            continue
                        if hap_name in hap_class:
                            hap_count = len(hap_dict[hap_name])
                            hap_name = hap_point_list[num][:length+1]
                            continue
                        if len(list(filterfalse(lambda x: not x.startswith(hap_name), hap_point_list))) < hap_count and len(list(filterfalse(lambda x: not x.startswith(hap_name), hap_point_list))) > 1:
                            hap_dict[hap_name] = [i for i in hap_point_list if i.startswith(hap_name)]
                            hap_class.append(hap_name)
                            hap_count = len(hap_dict[hap_name])
                            hap_name = hap_point_list[num][:length+1]
                        else:
                            hap_name = hap_point_list[num][:length+1]
            if hap_dict:
                level_list = [i]
                key_list = sorted(list(hap_dict.keys()))
                key_list.reverse()
                ele_list = []
                for key in key_list:
                    one_key_list = hap_dict[key]
                    ele_list.extend(one_key_list)
                ele_list.extend(key_list)
                ele_list = sorted(list(set(ele_list)))
                level_count = 0
                level_hap = i.split('\t')[-1]
                level_dict = {}
                level_dict[i.split('\t')[-1]] = 0
                for ele_num, ele in enumerate(ele_list[1:]):
                    if ele.startswith(level_hap):
                        level_count += 1
                        level_dict[ele] = level_count
                    else:
                        for ele_back in reversed(ele_list[:ele_num]):
                            if ele.startswith(ele_back):
                                level_count = level_dict[ele_back]+1
                                level_dict[ele] = level_count
                                break
                    level_list.append('-\t' * (level+level_count) + ele)
                    level_hap = ele
                level_list = sorted(level_list, key=lambda level_list: level_list.split('-\t')[-1])
                if i == level_list[0]:
                    new_tree_info.extend(level_list)
            else:
                new_tree_info.append(i)

        # remove tree trunk not for tree construction
        intersection = [i for i in new_tree_info if i.split('\t')[-1] not in filter_trunk]
        for i in reversed(intersection):
            if i.split('\t')[-1] in self.common_trunk:
                common_level = i.count('\t')
                if intersection.index(i) == len(intersection)-1:
                    filter_trunk.append(i.split('\t')[-1])
                elif intersection[intersection.index(i)+1].count('\t') != common_level+1:
                    filter_trunk.append(i.split('\t')[-1])
                    if i.split('\t')[-1] == 'P1':
                        filter_trunk.append('P')
                    elif i.split('\t')[-1] == 'NO1':
                        filter_trunk.append('NO')
                elif all([hap.split('\t')[-1] in self.common_trunk or len(hap.split('\t')[-1]) == 1 for hap in intersection[intersection.index(i):]]):
                    filter_trunk.extend([hap.split('\t')[-1] for hap in intersection[intersection.index(i):]])

        # remove duplicated branch in tree construction
        rm_list = []
        iter_num = 0
        for i in new_tree_info:
            node_hap = i.split('\t')[-1]
            node_level = i.count('\t')
            if node_hap in self.common_trunk or len(node_hap) == 1:
                continue
            if new_tree_info.index(i) == len(new_tree_info)-1:
                continue
            if node_hap in hap_df['Haplogroup'].tolist():
                continue
            if new_tree_info[new_tree_info.index(i)+1].split('\t')[-1].startswith(node_hap):
                branch_hap = new_tree_info[new_tree_info.index(i)+1].split('\t')[-1]
                node_base = [i for i in new_tree_info if i.split('\t')[-1].startswith(node_hap) and i.split('\t')[-1]!= node_hap]
                branch_base = [i for i in new_tree_info if i.split('\t')[-1].startswith(branch_hap)]
                if node_base == branch_base:
                    if iter_num == 0:
                        new_tree_info_irredundant = [(i.count('\t')-1)*'-\t' + i.split('\t')[-1] if i.split('\t')[-1].startswith(node_hap) else i for i in new_tree_info]
                    else:
                        new_tree_info_irredundant = [(i.count('\t')-1)*'-\t' + i.split('\t')[-1] if i.split('\t')[-1].startswith(node_hap) else i for i in new_tree_info_irredundant]
                    iter_num += 1
                    rm_list.append(node_hap)
        new_tree_info = [i for i in new_tree_info_irredundant if i.split('\t')[-1] not in rm_list]

        # add SampleID to tree
        ind_tree_info = []
        min_level_num = config.getint('PlotParameters', 'MinLevelNum')
        for i in new_tree_info:
            ind_hap = i.split('-\t')[-1]
            level_num = i.count('-\t')
            if ind_hap in hap_phylo['Haplogroup'].tolist():
                ind_tree_info.append(i)
                if isinstance(hap_inds.at[ind_hap, 'SampleID'], str):
                    ind_tree_info.append((level_num+1) * '-\t' + hap_inds.at[ind_hap, 'SampleID'])
                else:
                    ind_list = hap_inds.at[ind_hap, 'SampleID'].tolist()
                    ind_tree_info = ind_tree_info + [(level_num+1) * '-\t' + i for i in ind_list]
                if level_num < min_level_num:
                    min_level_num = level_num
            elif ind_hap not in filter_trunk:
                ind_tree_info.append(i)
        max_length = max([len(i.split('-\t')[-1]) for i in new_tree_info])

        # tree in newick format for visulization
        samples = hap_df['SampleID'].tolist()
        base, trunk_length = self._write_nwk(ind_tree_info, hap_df['SampleID'].tolist())
        base = base[1:].replace(');', ';') # root

        # remove reduandant nodes
        base_tree = Tree(base, format=8)
        for node in base_tree.traverse('preorder'):
            if (not node.is_leaf()) and len(node.children) == 1:
                child_name = node.children[0].name
                if child_name in self.common_trunk or len(child_name) == 1:
                    node.delete()

        base = base_tree.write(format=8)
        base = base.replace('NoName', '')

        return hap_df, base, trunk_length, max_length, samples

    # output tree in specific format
    def output_tree(self, base, tree_format):

        if base.endswith(');'):
            base = base[1:].replace(');', ';') # root

        # output tree in specific format
        output_file = self.output_set[1]
        if tree_format == 'newick':
            output_tree = open(output_file, 'w')
            output_tree.write(base)
            output_tree.close()
        else:
            from io import StringIO
            tmp_tree = StringIO(base)
            nwk_tree = Phylo.parse(tmp_tree, 'newick')
            Phylo.write(nwk_tree, output_file, tree_format)

    # apply sequence alignment method to construct a bifurcating tree
    def phylo_internal(self, ids, seqs, base, samples, align):

        # change the internal node of tree
        base_tree = Tree(base, format=8)
        for node in base_tree.iter_descendants('postorder'):
            if len(node.children) > 2 and node.name != '':
                node_name = node.name
                tmp_node = Tree('(tmp);', format=8)
                children = node.children
                children_haps_num = len([i for i in children if i.name not in samples])
                children_sample_num = len([i for i in children if i.name in samples])
                if children_haps_num == 0: # all children samples
                    children_samples = [i for i in children if i.name in samples]
                    self._add_single_tree_from_tmp(node, tmp_node, children_samples, node_name, ids, seqs, align, True)
                elif children_haps_num == 1:
                    if children_sample_num == 2:
                        children_samples_name = [i.name for i in children if i.name in samples]
                        children_haps_name = [i.name for i in children if i.name not in samples]
                        polynode = base_tree.get_common_ancestor(children_haps_name[0], children_samples_name[0])
                        polynode.resolve_polytomy(recursive=False)
                    else:
                        children_samples = [i for i in children if i.name in samples]
                        self._add_single_tree_from_tmp(node, tmp_node, children_samples, node_name, ids, seqs, align)
                elif children_haps_num > 1:
                    hap_sample_dict = {}
                    if children_sample_num <= 1: # all sub-haplogroups
                        children_haps = [i for i in children if i.name not in samples]
                        self._add_multi_tree_from_tmp(node, tmp_node, children_haps, node_name, ids, seqs, samples, align)
                    else:
                        children_samples = [i for i in children if i.name in samples]
                        self._add_single_tree_from_tmp(node, tmp_node, children_samples, node_name, ids, seqs, align)
                        children_haps = [i for i in children if i.name not in samples]
                        tmp_node2 = Tree('(tmp2);', format=8)
                        self._add_multi_tree_from_tmp(node, tmp_node2, children_haps, node_name, ids, seqs, samples, align)
            gc.collect()

        # remove nodes without children
        for node in base_tree.traverse('preorder'):
            if (not node.is_leaf()) and len(node.children) == 1:
                child_name = node.children[0].name
                if child_name not in samples:
                    node.delete()

        base = base_tree.write(format=8)
        base = base.replace('NoName', '')

        return base

    def _add_single_tree_from_tmp(self, node, tmp_node, children_samples, node_name, ids, seqs, align, replace=False):

        node.add_child(tmp_node)
        children_samples_name = [i.name for i in children_samples]
        # delete parallel samples
        for i in children_samples:
            i.delete()
        # concat sequence_distance_tree to haplogroup node
        id_idx = [ids.index(i) for i in children_samples_name]
        children_seqs = [seqs[i] for i in id_idx]
        children_tree = self._get_sequence_distance_tree(children_seqs, children_samples_name, align)
        if replace:
            children_tree = children_tree.replace(';', node_name+';')
        node.search_nodes(name='tmp')[0].add_child(Tree(children_tree, format=8))
        node.search_nodes(name='tmp')[0].delete()

    def _add_multi_tree_from_tmp(self, node, tmp_node, children_haps, node_name, ids, seqs, samples, align, replace=False):

        node.add_child(tmp_node)
        tmp_name = [i.name for i in tmp_node][0]
        hap_sample_dict = {}
        for hap in children_haps:
            hap_children_sample = list(filter(lambda x: x != '' and x in samples,
                                              [i.name for i in hap.traverse()]))[0]
            hap_sample_dict[hap_children_sample] = hap.name
        hap_children_samples = list(hap_sample_dict.keys())
        id_idx = [ids.index(i) for i in hap_children_samples]
        hap_children_seqs = [seqs[i] for i in id_idx]
        hap_children_tree = self._get_sequence_distance_tree(hap_children_seqs, hap_children_samples, align)
        for hap_children_sample, hap_name in hap_sample_dict.items():
            hap_children_tree = hap_children_tree.replace(hap_children_sample, hap_name+tmp_name)
        tmp_tree = Tree(hap_children_tree, format=8)
        for hap in children_haps:
            node_to_add = Tree(hap.write(format=8), format=8)
            tmp_tree.search_nodes(name=hap.name+tmp_name)[0].add_child(node_to_add)
            tmp_tree.search_nodes(name=hap.name+tmp_name)[0].delete()
            hap.detach()
        node.search_nodes(name=tmp_name)[0].add_child(tmp_tree)
        node.search_nodes(name=tmp_name)[0].delete()

    # get sub-tree from multiple sequence alignment
    def _get_sequence_distance_tree(self, seqs, ids, align):

        new_seqs = []
        # prune sequences
        for i in range(len(seqs[0])):
            all_al = [s[i] for s in seqs]
            al = list(set(all_al))
            N_count = 0
            if 'N' in al:
                al.remove('N')
                N_count += all_al.count('N')
            if '-' in al:
                al.remove('-')
                N_count += all_al.count('-')
            N_rate = N_count / len(all_al)
            if len(al) != 1 and N_rate < 0.05:
                new_seqs.append(all_al)

        new_seqs = np.array(new_seqs).T.tolist()

        # IBS method
        if align == 'ibs':
            from scipy.cluster import hierarchy
            from scipy.spatial.distance import squareform

            pairwise_df = pd.DataFrame(index=ids, columns=ids)

            # get pairwise IBS
            ind_num = 0
            for ind1 in ids:
                for ind2 in ids[ind_num:]:
                    if ind1 == ind2:
                        pairwise_df.at[ind1, ind2] = 1
                    else:
                        ibs = self._calculate_ibs(seqs, ids, ind1, ind2)
                        pairwise_df.at[ind1, ind2] = ibs
                        pairwise_df.at[ind2, ind1] = ibs
                ind_num += 1
            pairwise_df = 1 - pairwise_df

            # convert distance matrix to hierarchy
            dm = hierarchy.linkage(squareform(np.array(pairwise_df)), 'average')
            dm_tree = hierarchy.to_tree(dm, False)
            children_tree = self._hier_to_newick(dm_tree, '', pairwise_df.index.tolist())
        else:
            import tempfile
            from Bio import Phylo
            from Bio.Seq import Seq
            from Bio.SeqRecord import SeqRecord
            from Bio.Align import MultipleSeqAlignment
            from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

            new_seqs = [''.join(i) for i in new_seqs]
            seq_list = [SeqRecord(Seq(j), id=i, name=i) for i, j in zip(ids, new_seqs)]
            aln = MultipleSeqAlignment(seq_list)
            calculator = DistanceCalculator('identity')

            # UPGMA method
            if align == 'upgma':

                constructor = DistanceTreeConstructor(calculator, 'upgma')

                tree = constructor.build_tree(aln)
                tree_tmp_file = tempfile.NamedTemporaryFile(mode='w')
                Phylo.write(tree, tree_tmp_file.name, 'newick')
                ete_tree = Tree(open(tree_tmp_file.name).read(), format=1)
                children_tree = ete_tree.write(format=9)

            # MP (maximum parsimony) method
            elif align == 'mp':

                from Bio.Phylo.TreeConstruction import ParsimonyScorer, NNITreeSearcher, ParsimonyTreeConstructor

                dm = calculator.get_distance(aln)
                constructor = DistanceTreeConstructor(calculator, 'nj')
                tree = constructor.build_tree(aln)
                scorer = ParsimonyScorer()
                searcher = NNITreeSearcher(scorer)
                constructor = ParsimonyTreeConstructor(searcher, tree)
                pars_tree = constructor.build_tree(aln)

                tree_tmp_file = tempfile.NamedTemporaryFile(mode='w')
                Phylo.write(pars_tree, tree_tmp_file.name, 'newick')
                ete_tree = Tree(open(tree_tmp_file.name).read(), format=1)
                children_tree = ete_tree.write(format=9)

        return children_tree

    # IBS calculation
    def _calculate_ibs(self, seqs, ids, sample1, sample2):

        index1 = ids.index(sample1)
        index2 = ids.index(sample2)
        ibs0 = 0
        ibs1 = 0
        ibs2 = 0
        for i, j in zip(seqs[index1], seqs[index2]):
            if i in ['N', '-'] or j in ['N', '-']:
                continue
            if i == j:
                ibs2 += 1
            else:
                ibs0 += 1
        ibs = (ibs2) / (ibs0+ibs2)

        return ibs

    # convert hierarchy to newick
    def _hier_to_newick(self, node, newick, leaf_names):

        if node.is_leaf():
            return '%s%s' % (leaf_names[node.id], newick)
        else:
            if len(newick) > 0:
                newick = ')%s' % (newick)
            else:
                newick = ');'
            newick = self._hier_to_newick(node.get_left(), newick, leaf_names)
            newick = self._hier_to_newick(node.get_right(), ',%s' % (newick), leaf_names)
            newick = '(%s' % (newick)

        return newick

    # function used for construction of write newick file
    def _write_nwk(self, ind_tree_info, samples):

        base = ''
        width, level, length =0, 0, 0.5
        trunk_list = []
        trunk_length = {}
        count_num = ind_tree_info[-1].count('-')
        for i in ind_tree_info[1:]:
            trunk = i.split('\t')[-1]
            width_new = i.count('-\t')
            if width_new > width:
                width = width_new
            level_new = i.count('-')-trunk.count('-')
            trunk_list.append(trunk)
            dif = level_new - level
            if dif == 0:
                length = length
                base = trunk + ',' + base
                if ind_tree_info.index(i)+1 == len(ind_tree_info):
                    base =  '('*(base.count(')')-base.count('(')) + base
            elif dif > 0:
                length = 1+1*(count_num-level_new)
                base = trunk + ')' + base
                if ind_tree_info.index(i)+1 == len(ind_tree_info):
                    base =  '('*(base.count(')')-base.count('(')) + base
            elif dif < 0:
                length = 1+1*(count_num-level_new)
                if ind_tree_info.index(i)+1 != len(ind_tree_info):
                    base =  trunk + ',' + '('*abs(dif) + base
                else:
                    base =  '('*(base.count(')')-base.count('(')-1) + trunk + ',(' + base
            trunk_length[trunk] = length
            level = level_new

        base = '(' + base + 'Y-Adam);'

        return base, trunk_length

    # plot tree figure
    def _visualize_tree(self, base, trunk_length, max_length, samples, layout, plot_info, color_dict=None):

        from ete3 import TreeStyle, NodeStyle, TextFace, RectFace
        # build tree from newick
        base_tree = Tree(base, format=8)

        num = 0
        for node in base_tree.traverse('postorder'):
            if node.name == '':
                node.name = 'BLANK_NODE%d' % num
            num += 1

        # add text in trunk and leaf, and keep all branch length same in topology
        for node in base_tree.traverse('postorder'):
            trunk = base_tree & node.name
            if node.name in samples:
                text = TextFace(node.name, ftype='Arial', fsize=18)
                text.margin_left = 4
                trunk.add_face(text, column=0, position='branch-right')
            elif node.name.startswith('BLANK_NODE'):
                text = ''
                trunk.add_face(TextFace(text, ftype='Arial', fsize=18), column=0, position='branch-top')
            else:
                space_num = max_length - len(node.name)
                text = node.name + ' '*space_num
                trunk.add_face(TextFace(text, ftype='Arial', fsize=18), column=0, position='branch-top')
            if node.is_leaf():
                node.dist = trunk_length[node.name]
            else:
                node.dist = 1
            gc.collect()

        # parameters for node style
        # the only difference between TreeNode.traverse() and TreeNode.iter_descendants() is that the first will include the root node in the iteration
        node_color = {}
        if isinstance(plot_info, pd.DataFrame): # match colors to population
            label = plot_info.columns[0]
            for node in base_tree.traverse('postorder'):
                ns = NodeStyle()
                ns['hz_line_width'] = 1
                ns['vt_line_width'] = 1
                if node.is_leaf() and node.name in samples:
                    # color
                    color = color_dict[plot_info.at[node.name, label]]
                    ns['fgcolor'] = color
                    node_color[node.name] = color
                    # shape
                    ns['shape'] = 'square'
                    ns['size'] = 50
                else:
                    children_nodes = [c_node for c_node in node.children]
                    node_colors = [node_color[i.name] for i in children_nodes]
                    if len(set(node_colors)) > 1:
                        if node.name in self.common_trunk or len(node.name) == 1:
                            color = 'black'
                        else:
                            for g_node in children_nodes:
                                node_colors += [node_color[i.name] for i in g_node]
                    node_color[node.name] = color
                    ns['fgcolor'] = 'black'
                    ns['size'] = 0
                ns['vt_line_color'] = color
                ns['hz_line_color'] = color
                node.set_style(ns)
        else: # all the leaf is same color
            for node in base_tree.traverse():
                ns = NodeStyle()
                ns['hz_line_width'] = 1
                ns['vt_line_width'] = 1
                ns['fgcolor'] = 'black'
                ns['size'] = 0
                node.set_style(ns)

        # parameters for tree style
        ts = TreeStyle()
        ts.min_leaf_separation = 30
        ts.show_leaf_name = False
        ts.show_branch_length = False
        ts.scale = max_length*15
        ts.branch_vertical_margin = 0
        ts.optimal_scale_level = 'full'
        ts.show_scale = False
        ts.margin_right = 20
        ts.margin_left = 20
        if layout == 'circular':
            ts.mode = 'c'

        # add legend and output tree
        base_tree.convert_to_ultrametric()
        if color_dict:
            if label == 'Color':
                # output figure with color information only
                base_tree.render(self.output_set[0]+'.pdf', tree_style=ts)
            else:
                # add legend
                for name, color in color_dict.items():
                    icon_legend = RectFace(width=200, height=200, fgcolor='black', bgcolor=color)
                    icon_legend.margin_left = 200
                    icon_legend.margin_right = 100
                    icon_legend.margin_top = 50
                    icon_legend.margin_bottom = 50
                    ts.legend.add_face(icon_legend, column=0)
                    text_legend = TextFace(name, fsize=200, ftype='Arial')
                    text_legend.margin_top = 50
                    text_legend.margin_bottom = 50
                    ts.legend.add_face(text_legend, column=1)
                ts.legend_position = 1
                # output figure with population information
                base_tree.render(self.output_set[0]+'.%s.pdf' % label, tree_style=ts)
        else:
            # output figure without population
            base_tree.render(self.output_set[0]+'.pdf', tree_style=ts)

    # match color to populations
    def _match_colors(self, population_info, label):

        # create a empty dict to match color and population
        color_dict = {}
        # if preset color, check color names
        if label == 'Color':
            colors = population_info[label]
            for i in set(colors):
                if matplotlib.colors.is_color_like(i):
                    color_dict[i] = i
                else:
                    self.logger.error('[Y-LineageTracker] [Error] %s is not a RGB color' % i)
                    sys.exit()
        else:
            # match color to populations
            populations = sorted(list(set(population_info[label])))
            from FilesIO import set_color_num
            color_num = len(populations)
            colors, cmap = set_color_num(color_num, get_cmap=True)

            # creat a dict usef for match color to population
            for i, j in zip(populations, cmap):
                color_dict[i] = matplotlib.colors.to_hex(j)

        plot_info = population_info[[label]]

        return plot_info, color_dict

    # main function for tree plot
    def output_tree_figure(self, population_info, base, trunk_length, max_length, samples, layout):

        if isinstance(population_info, pd.DataFrame):
            if 'Color' in population_info.columns:
                plot_info, color_dict = self._match_colors(population_info, 'Color')
                self._visualize_tree(base, trunk_length, max_length, samples, layout, plot_info, color_dict)
            else:
                if 'Population' in population_info.columns:
                    plot_info, color_dict = self._match_colors(population_info, 'Population')
                    self._visualize_tree(base, trunk_length, max_length, samples, layout, plot_info, color_dict)
                if 'Group' in population_info.columns:
                    plot_info, color_dict = self._match_colors(population_info, 'Group')
                    self._visualize_tree(base, trunk_length, max_length, samples, layout, plot_info, color_dict)
        else:
            plot_info = None
            self._visualize_tree(base, trunk_length, max_length, samples, plot_info, layout)


def main():

    start = time.perf_counter()
    arguments = phylo_parser()

    # set of output files
    output_set = check_phylo_output(arguments.tree, arguments.output)

    from FilesIO import check_hap_input, check_population, check_overlap, time_count
    # check haplogroup data
    hap_data = check_hap_input(arguments.hg, 'haplogroup')

    # check population data
    if arguments.population:
        population_data = check_population(arguments.population)
        hap_data, population_data = check_overlap(hap_data, population_data)

    # set log file
    set_log(arguments, output_set[-1])

    # start phylogeny analysis
    phylo = PhyloHaplogroup(hap_data, output_set)
    # read sequence alignment file
    if arguments.seq:
        ids, seqs = phylo.read_seq_file(arguments.seq, arguments.format)
        hap_df, base, trunk_length, max_length, samples = phylo.get_base_tree(ids) # build base tree
    else:
        hap_df, base, trunk_length, max_length, samples = phylo.get_base_tree()
    # construct sub-tree internal nodes
    if arguments.seq:
        base = phylo.phylo_internal(ids, seqs, base, samples, arguments.align)
    # combine with population data
    if arguments.population:
        population_info = pd.merge(hap_df, population_data, on='SampleID')
        population_info = population_info.set_index('SampleID')
    else:
        population_info = None
    # output tree in specific format
    phylo.output_tree(base, arguments.tree)
    # output tree figure
    phylo.output_tree_figure(population_info, base, trunk_length, max_length, samples, arguments.layout)

    time_count(start)


if __name__ == '__main__':
    main()
