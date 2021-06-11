import os
import sys
import time
import logging
import argparse
import pandas as pd

sys.path.insert(0, '..')
from GetLineage import MatchVariant, TrackLineage
from FilesIO import CommonData

common_data = CommonData()
test_path = os.path.split(os.path.realpath(__file__))[0] + '/TestData/'
out_path = test_path + 'TestOut/'
if not os.path.exists(out_path):
    os.makedirs(out_path)


def test_classify(input, format, build, sample):

    cutoff = None
    snp_only = False
    appro = False
    chip = False
    mutation = True
    output_set = [out_path + 'Test.hapresult.hg',
                  out_path + 'Test.lineageresult.hg',
                  out_path + 'Test.mutationresult.txt']

    # read ISOGG haplogroup reference data
    ref_info, ref_ind_info = common_data.read_ref_info(snp_only, appro)
    # filter samples
    from ProcessData import ReadFilterData
    read_filter = ReadFilterData.ReadFilterData(input, sample)
    if format == 'bam':
        data_info, inds = read_filter.read_bam(ref_info, build, 'snp')
    else:
        data_info, inds = read_filter.read_filter_variant(cutoff, format)
    # match input data to ISOGG reference
    match = MatchVariant.MatchVariant(build)
    matched_info, header_cut = match.match_variant_snp(ref_info, data_info, format)

    # classify haplogroups and build lineage track
    tracker = TrackLineage.TrackLineage(matched_info, ref_ind_info, build, format, chip)
    hap_data = tracker.classify_haplogroup(inds, matched_info, ref_info, output_set, mutation, header_cut)


def test_cluster(hg_data, pop_data, method):

    from ClusterHaplogroup import ClusterAnalysis

    # set parameters
    level = 3
    freq = True
    if method == 'pca':
        output_set = [out_path + 'Test.PCA.fig.pdf',
                      out_path + 'Test.PCA.eigenvectors.txt',
                      out_path + 'Test.PCA.eigenval.txt',
                      out_path + 'Test.level.hg',
                      out_path + 'Test.freq.txt']
    else:
        output_set = [out_path + 'Test.MDS.fig.pdf',
                      out_path + 'Test.MDS.fst.txt',
                      out_path + 'Test.MDS.embedding.txt',
                      out_path + 'Test.level.hg',
                      out_path + 'Test.freq.txt']

    # start cluserting analysis
    cluster = ClusterAnalysis()
    hg_data = cluster.classify_hap_data(level, hg_data, output_set[3])
    cluster.output_cluster(method, freq, output_set, hg_data, pop_data)


def test_phylo(hg_data, pop_data, seq_file):

    from PhyloHaplogroup import PhyloHaplogroup

    align = 'mp'
    format = 'vcf'
    tree = 'newick'
    layout = 'rectangular'

    # set of output files
    output_set = [out_path + 'Test.phylo.fig.pdf',
                  out_path + 'Test.tree.nwk']

    # phylo analysis
    phylo = PhyloHaplogroup(hg_data, output_set)
    ids, seqs = phylo.read_seq_file(seq_file, format)
    hap_df, base, trunk_length, max_length, samples = phylo.get_base_tree(ids)
    base = phylo.phylo_internal(ids, seqs, base, samples, align)
    pop_info = pd.merge(hap_df, pop_data, on='SampleID')
    pop_info = pop_info.set_index('SampleID')
    phylo.output_tree(base, tree)
    #phylo.output_tree_figure(pop_info, base, trunk_length, max_length, samples, layout)


def test_genostr():

    import pysam


def test_net(hap_data, pop_data):

    from HaploNetwork import HaploNetwork

    class NetArgs(object):

        gap = 'allele'
        filter = 0.7
        treetype = 'mjn'
        search = 0
        alt = 'hide'
        mediansize = 0
        fdi = True
        euqalsize = False

    arguments = NetArgs()
    gap = 'allele'
    label = 'Population'
    path = out_path + 'Test'

    # calculate network and visualize result
    net = HaploNetwork(gap, path)
    net.calculate_and_output(hap_data, pop_data, arguments, label)


def test_stat(hap_data, pop_data):

    from StatisticsAnalysis import StatisticsAnalysis, predict_haplogroup

    class StatArgs(object):

        hd = True
        mpd = True
        gd = False
        amova = 'fst'

    arguments = StatArgs()
    path = out_path + 'Test'
    label = 'Population'
    stat = StatisticsAnalysis(label, hap_data, pop_data)
    stat.calculate_and_output(arguments, path)
    predict_haplogroup(hap_data, path)


def test_time(seq_file, tree_file):

    from EstimateTime import EstimateTime

    class TimeArgs(object):

        tree = tree_file
        seq = seq_file
        format = 'vcf'
        filter = 0.1
        mcmc = 2000
        rate = 7.6e-10
        model = 'hky85'
        calibration = 'O2a:10000-13000'
        auto = False

    arguments = TimeArgs()
    path = out_path + 'Test'

    # estimate divergence time
    estimate = EstimateTime(path)
    estimate.run_mcmc(arguments)

def test_tmrca(input_file, format, pop_file):

    from CommonAncestor import CommonAncestor

    class TmrcaArgs(object):

        if format == 'net':
            net = input_file
            matrix = None
            population = None
        else:
            net = None
            matrix = input_file
            population = pop_file
        rate = 6e-3
        generation = 30
        ad = ['NA18530:NA18543,NA18544']

    arguments = TmrcaArgs()
    path = out_path + 'Test'

    estimate = CommonAncestor()
    estimate.run_tmrca(arguments, path)


def main():

    from FilesIO import check_hap_input, check_population, check_overlap, set_haplotype_index

    vcf_file = test_path + '/KGP_sample7.vcf'
    tree_file = test_path + '/KGP_sample7.tree.nwk'
    net_file = test_path + '/KGP_sample7.net'

    hg_file = test_path + '/KGP_hgs.hg'
    hg_data = check_hap_input(hg_file, 'haplogroup')

    pop_file = test_path + '/KGP_pops.pop'
    pop_data = check_population(pop_file)

    hap_file = test_path + '/KGP_O2_sample7.matrix'
    hap_data = check_hap_input(hap_file, 'haplotype', 'matrix')

    hap_data, pop_data2 = check_overlap(hap_data, pop_data)
    hap_data = set_haplotype_index(hap_data, 'matrix')

    print('[Y-LineageTracker] Testing classify command...')
    try:
        test_classify(vcf_file, 'vcf', 37, None)
        print('[Y-LineageTracker] Successfully tested classify command')
    except:
        print('[Y-LineageTracker] Failed in testing classify command')

    print('[Y-LineageTracker] Testing cluster command...')
    try:
        test_cluster(hg_data, pop_data, 'pca')
        test_cluster(hg_data, pop_data, 'mds')
        print('[Y-LineageTracker] Successfully tested clster command')
    except:
        print('[Y-LineageTracker] Failed in testing cluster command')

    print('[Y-LineageTracker] Testing phylo command...')
    try:
        test_phylo(hg_data, pop_data, vcf_file)
        print('[Y-LineageTracker] Successfully tested phylo command')
    except:
        print('[Y-LineageTracker] Failed in testing phylo command')

    print('[Y-LineageTracker] Testing genostr command...')
    try:
        test_genostr()
        print('[Y-LineageTracker] Successfully tested genostr command')
    except:
        print('[Y-LineageTracker] Failed in testing genostr command')

    print('[Y-LineageTracker] Testing net command...')
    try:
        test_net(hap_data, pop_data2)
        print('[Y-LineageTracker] Successfully tested net command')
    except:
        print('[Y-LineageTracker] Failed in testing net command')

    print('[Y-LineageTracker] Testing stat command...')
    try:
        test_stat(hap_data, pop_data2)
        print('[Y-LineageTracker] Successfully tested stat command')
    except:
        print('[Y-LineageTracker] Failed in testing stat command')

    print('[Y-LineageTracker] Testing time command...')
    try:
        test_time(vcf_file, tree_file)
        print('[Y-LineageTracker] Successfully tested time command')
    except:
        print('[Y-LineageTracker] Failed in testing time command')

    print('[Y-LineageTracker] Testing tmrca command...')
    try:
        test_tmrca(hap_file, 'hap', pop_file)
        test_tmrca(net_file, 'net', pop_file)
        print('[Y-LineageTracker] Successfully tested tmrca command')
    except:
        print('[Y-LineageTracker] Failed in testing tmrca command')


if __name__ == '__main__':
    main()
