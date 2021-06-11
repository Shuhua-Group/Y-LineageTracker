import sys
import gzip
import logging
import tempfile
import itertools
import numpy as np
import pandas as pd

from GetConfig import getConfig
from collections import Counter


config = getConfig()


class ReadFilterData(object):
    '''
    This class is mainly used to:
    1. read VCF and BAM file
    2. get individuals in sample list
    3. filter female individual by missing rate or K-means
    4. filter sites to get region for time estimation
    5. genotype Y-STRs after read BAMs
    6. map ISOGG reference to BAMs, and get a list for subsequent analysis
    '''

    def __init__(self, input_file, samples):

        self.logger = logging.getLogger()
        self.samples = samples
        self.input_file = input_file

    # get overlaped male individuals from sample list
    def _get_male_inds(self, inds):

        males_not_in_data = [i for i in inds if i not in self.samples]
        if males_not_in_data:
            if males_not_in_data == self.samples:
                self.logger.error('[Y-LineageTracker] [Error] No samples in input data, please check sample file')
                used_males = inds
            else:
                self.logger.warning('[Y-LineageTracker] [Warning] %s not in input data' % (', '.join(males_not_in_data)))
                used_males = [i for i in inds if i in self.samples]
        else:
            used_males = self.samples

        return used_males

    # get motif number of Y-STR
    def _count_motif_num(self, bam, contig, start, end, base_motif_num, motif, motif_length, end_motif):

        add_end = end+motif_length
        pos_info = bam.fetch(contig, start-1, add_end)
        add_ref_length = add_end-start+1

        # get all possible sequence of specific region
        seqs = []
        for i in pos_info:
            num = start-i.pos-1
            num2 = add_end-i.pos
            if num > 0:
                if num2 < i.qlen:
                    seq = i.seq[num:num2]
                else:
                    seq_tmp = i.seq[num: i.qlen]
                    seq = seq_tmp+(add_ref_length-len(seq_tmp))*'.'
            else:
                seq_tmp = i.seq[:num2]
                seq = (add_ref_length-len(seq_tmp))*'.'+seq_tmp
            seqs.append(seq)

        # convert seq region to block
        block = [''.join(s) for s in zip(*seqs)]
        # get final sequence from block
        final_seq = ''.join([Counter(i.strip('.')).most_common()[0][0] if i.count('.') != len(i) else '.' for i in block])

        if len(final_seq) == 0 or '.' in final_seq:
            sample_motif_count = '.'
        else:
            if len(motif) == 1:
                sample_motif_count = final_seq.count(motif[0])
            else:
                sample_motif_count = sum([final_seq.count(m) for m in motif])

            # recursion if motif+1
            if sample_motif_count > base_motif_num:
                if pd.isna(end_motif) or (not pd.isna(end_motif) and final_seq.endswith(end_motif)):
                    base_motif_num = sample_motif_count
                    end += motif_length
                    sample_motif_count = self._count_motif_num(bam, contig, start, end, base_motif_num, motif, motif_length, end_motif)

        return sample_motif_count

    def _get_specific_sites(self, build, sub_ref_info, id_info_dict):


        data_info_array = []
        for n, pos_num in enumerate(sub_ref_info['Build'+str(build)]):
            line = []
            int_pos = int(pos_num)
            for bam, contig in id_info_dict.values():
                pos_info = bam.fetch(contig, int_pos-1, int_pos)
                pos_list = [i.seq[int_pos-i.pos-1].upper() for i in pos_info if int_pos-i.pos < i.qlen]
                if len(set(pos_list)) == 1:
                    pos = pos_list[0]
                elif len(set(pos_list)) > 1:
                    pos = Counter(pos_list).most_common()[0][0]
                else:
                    pos = '.'
                line.append(pos)
            data_info_array.append([pos_num] + line)

        data_info = pd.DataFrame(columns=['POS']+list(id_info_dict.keys()), data=np.array(data_info_array))

        return data_info

    def _match_sites(self, data_info, ref_info, build, matched_marker, header_cut, header_end):

        matched_info = pd.merge(data_info, ref_info,
                                left_on='POS',
                                right_on='Build'+str(build),
                                how='inner',
                                sort=False)
        matched_info_array = np.array(matched_info)
        hap_function = (lambda x: matched_marker+matched_info_array[i, -10] if x == match_marker else x)
        for i in range(len(matched_info)):
            match_marker = matched_info_array[i, -6]
            matched_info_array[i, header_cut: header_end] = list(map(hap_function, matched_info_array[i, header_cut: header_end]))

        final_matched_info = pd.DataFrame(columns=matched_info.columns.tolist(), data=matched_info_array)

        return final_matched_info

    def _get_female_by_KMeans(self, missing_rate_data, samples):

        from sklearn.cluster import KMeans

        missing_rate_data.extend([0, 1])
        missing_rate_data = np.array(missing_rate_data).reshape(-1, 1)
        estimator = KMeans(n_clusters=2)
        res = estimator.fit_predict(missing_rate_data)[: -2]
        centroids = estimator.cluster_centers_
        if centroids[0][0] > centroids[1][0]:
            del_list = [samples[j] for i, j in zip(res, range(len(res))) if i == 0]
            female_num = len(list(filter(lambda x: x == 0, res)))
        else:
            del_list = [samples[j] for i, j in zip(res, range(len(res))) if i == 1]
            female_num = len(list(filter(lambda x: x == 1, res)))

        return del_list, female_num

    # read BAM and do matching analysis
    def read_bam(self, ref_info, build, extract_type, format):

        import os
        import pysam
        from collections import Counter

        if format == 'bam':
            open_method = 'rb'
        else:
            open_method = 'rc'

        # get bam list
        input_file_list = self.input_file.split(', ')
        if len(input_file_list) == 1:
            file = input_file_list[0]
            try:
                pysam.AlignmentFile(file, open_method)
                bam_list = input_file_list
            except ValueError:
                bam_list = open(self.input_file).read().splitlines()
        else:
            bam_list = self.input_file.split(', ')
            for i in bam_list:
                file_exists = os.path.isfile(i)
                if not file_exists:
                    raise FileNotFoundError('No such file or directory: %s' % i)

        id_info_dict = {}
        format_contigs = config.get('DataInfo', 'Contigs').split(',')
        for f in bam_list:
            file_name = os.path.basename(f)
            bam = pysam.AlignmentFile(f, open_method)
            bam_with_index = bam.has_index()
            if not bam_with_index:
                self.logger.error('[Y-LineageTracker] bam/cram file %s has no index' % f)
                sys.exit()
            ind = bam.header.get('RG')[0]['ID']
            if ind in self.samples:
                continue
            sample_contigs = [i.contig for i in bam.get_index_statistics()]
            overlaped_contigs = list(set(sample_contigs) & set(format_contigs))
            if len(overlaped_contigs) == 1:
                contig = overlaped_contigs[0]
            else:
                if len(overlaped_contigs) > 1:
                    self.logger.error('[Y-LineageTracker] More than one contigs for Y chromosome (Y, chrY, 24, chr24), please check input data')
                    sys.exit()
                if len(overlaped_contigs) == 0:
                    self.logger.error('[Y-LineageTracker] No contig for Y chromosome (Y, chrY, 24, chr24), please check input data')
                    sys.exit()
            id = ind+'('+file_name+')'
            id_info_dict[id] = [bam, contig]

        # for NRY haplogroup classification
        # match ISOGG sites to BAM files
        if extract_type == 'snp':
            # match sites of main trunk
            matched_marker = config.get('HaplogroupMarker', 'MatchMarker') # used to mark matched haplogroups
            header_cut = config.getint('DataInfo', 'BamHeaderCut')
            header_end = -config.getint('DataInfo', 'ISOGGInfoCut')
            common_trunk = config.get('HaplogroupTree', 'CommonTrunk').split(',')
            special_haps = config.get('HaplogroupTree', 'SpecialTrunk').split(',')

            # get allele of haplogroup on main trunks
            main_trunk_ref_info = ref_info[ref_info['MainInfo']=='Main']
            self.logger.info('[Y-LineageTracker] Extrating mian trunks from BAMs/CRAMs...')
            main_trunk_info = self._get_specific_sites(build, main_trunk_ref_info, id_info_dict)
            # filter females
            main_marker_num = main_trunk_info.index.size
            f_samples = main_trunk_info.columns.tolist()[1:]
            missing_rate_data = [main_trunk_info[col].tolist().count('.')/main_marker_num for col in f_samples]
            del_list, female_num = self._get_female_by_KMeans(missing_rate_data, f_samples)
            if del_list:
                if female_num == len(f_samples):
                    self.logger.error('[Y-LineageTracker] Program stopped since no male sample left for subsequent analysis')
                    sys.exit()
                else:
                    main_trunk_info.drop(del_list, inplace=True, axis=1)
                    for female in del_list:
                        del id_info_dict[female]

            # match main haplogroups
            main_matched_info = self._match_sites(main_trunk_info, main_trunk_ref_info, build, matched_marker, header_cut, header_end)
            self.logger.info('[Y-LineageTracker] Extrating mian trunks from BAMs/CRAMs finifhsed')

            # find which main haplogroup with higher matching rate
            main_sub_info = main_matched_info[main_matched_info['Haplogroup'].map(lambda x: x not in common_trunk)]
            all_haps = main_sub_info['Haplogroup'].tolist()
            main_sub_info = np.array(main_sub_info)
            id_hap_dict = {}
            for id_num, id in enumerate(id_info_dict.keys()):
                col_info = main_sub_info[:, id_num+header_cut]
                id_haps = sorted(list(map(lambda x: x[1:], list(filter(lambda x: x[0]==matched_marker, col_info)))))
                ratio = 0
                for hap in set(id_haps):
                    ratio2 = id_haps.count(hap) / all_haps.count(hap)
                    if ratio2 > ratio:
                        ratio = ratio2
                        main_hap = hap
                if ratio < 0.25:
                    main_hap = None
                id_hap_dict[id] = main_hap

            # extract sub-haplogroups of each main trunk
            self.logger.info('[Y-LineageTracker] Extrating sub-trunks from BAMs/CRAMs...')
            sub_trunk_ref_info = ref_info[ref_info['MainInfo']!='Main']
            info_list = []
            NA_ids = []
            for hap in sorted(set(id_hap_dict.values())):
                hap_ids = list({k: v for k, v in id_hap_dict.items() if v==hap}.keys())
                sub_id_info_dict = {k: v for k, v in id_info_dict.items() if k in hap_ids}
                if main_hap:
                    hap_ref_info = sub_trunk_ref_info[sub_trunk_ref_info['Haplogroup'].map(lambda x: x.startswith(hap))]
                    hap_info = self._get_specific_sites(build, hap_ref_info, sub_id_info_dict)
                    info_list.append(hap_info)
                else:
                    hap_ref_info = sub_trunk_ref_info[sub_trunk_ref_info['Haplogroup'].map(lambda x: x.startswith('A') and x in special_haps)]
                    hap_info = self._get_specific_sites(build, hap_ref_info, sub_id_info_dict)
            sub_hap_info = pd.concat(info_list, axis=0, sort=False).fillna('.')
            self.logger.info('[Y-LineageTracker] Extrating sub-trunks from BAMs/CRAMs finifhsed')

            sub_trunk_info = sub_hap_info[['POS']+list(id_info_dict.keys())]
            sub_matched_info = self._match_sites(sub_trunk_info, sub_trunk_ref_info, build, matched_marker, header_cut, header_end)

            data_info = np.array(pd.concat([main_matched_info, sub_matched_info]).drop_duplicates())
            inds = list(id_info_dict.keys())

            self.logger.info('[Y-LineageTracker] Extrating information from BAMs/CRAMs finifhsed')

        # for Y-STR genotyping
        # genotype Y-STRs from BAMs
        else:
            from FilesIO import get_str_from_panel
            str_info = get_str_from_panel(ref_info)
            data_info = pd.DataFrame(columns=list(id_info_dict.keys()))
            for num, STR_idx in enumerate(str_info.index):
                # read basic information
                start = str_info.at[STR_idx, 'Start'+str(build)]
                end = str_info.at[STR_idx, 'End'+str(build)]
                motif = str_info.at[STR_idx, 'RefMotif'].split(', ')
                str_name = str_info.at[STR_idx, 'Name'+str(build)]
                motif_length = str_info.at[STR_idx, 'MotifLength']
                end_motif = str_info.at[STR_idx, 'End']
                motif_count = []
                # matching
                for bam, contig in id_info_dict.values():
                    base_motif_num = str_info.at[STR_idx, 'Ref'+str(build)]
                    sample_motif_count = self._count_motif_num(bam, contig, start, end, base_motif_num, motif, motif_length, end_motif)
                    motif_count.append(sample_motif_count)
                data_info.loc[str_name] = motif_count
                percent = float(num + 1) * 100 / float(str_info.index.size)
                sys.stdout.write('[Y-LineageTracker] Genotyping Y-STR…… %.2f' % percent)
                sys.stdout.write('%\r')
                sys.stdout.flush()
            data_info = data_info.T
            inds = list(id_info_dict.keys())

            print()
            self.logger.info('[Y-LineageTracker] Genotyping Y-STR finished')

        return data_info, inds

    # open input file of variant calling
    def _read_variant(self, file_format):

        # open VCF
        if file_format == 'gzvcf':
            opened_file = gzip.open(self.input_file, 'rt')
        else:
            opened_file = open(self.input_file, 'r')

        # get row number of header
        config = getConfig()
        if  'vcf' in file_format:
            header_cut = config.getint('DataInfo', 'VcfHeaderCut')
            start_symbol = '#CHROM'
            chr_symbol = start_symbol
            for i in itertools.count(start=0, step=1):
                vcf_line = opened_file.readline()
                if vcf_line.startswith(start_symbol):
                    head_num = i
                    break
        else:
            header_cut = config.getint('DataInfo', 'InpHeaderCut')
            start_symbol = 'dbSNP'
            blank_num = 0
            for i in itertools.count(start=0, step=1):
                inp_line = opened_file.readline()
                if inp_line == '\n':
                    blank_num += 1
                if inp_line.startswith(start_symbol):
                    head_num = i - blank_num
                    start_symbol = inp_line.split('\t')[0]
                    chr_symbol = inp_line.split('\t')[1]
                    pos_symbol = inp_line.split('\t')[2]
                    break
        opened_file.close()

        # open input file in pandas
        self.logger.info('[Y-LineageTracker] Reading input file...')
        data = pd.read_csv(self.input_file,
                           header=head_num,
                           sep='\t',
                           encoding='utf-8',
                           dtype='object')
        self.logger.info('[Y-LineageTracker] Input file Read')

        # change header name of site vcf if input is INP
        if file_format == 'inp':
            data.rename(columns={pos_symbol: 'POS'}, inplace=True)
        # chr column keeps only Y/24
        chr_type = len(set(data[chr_symbol].tolist()))
        if chr_type > 1:
            print('[Y-LineageTracker] [Warning] Other chromosomes are detected in inputfile, only the Y chromosome is reserved for analysis')
            data = data[data[chr_symbol].map(lambda x: 'Y' in x or '24' in x or 'chrY' in x or 'chr24' in x)]

        return data, header_cut

    # filter individuals
    def _filter_variant(self, data, cutoff, file_format, header_cut):

        # list of indivisuals and number
        inds = data.columns.tolist()[header_cut:]
        ind_num_before = len(inds)
        self.logger.info('[Y-LineageTracker] There are %d individuals in input file' % (ind_num_before))

        # filter according to sample list
        if self.samples:
            used_males = self._get_male_inds(inds)
            header = data.columns.tolist()[:header_cut]
            data = data[header+used_males]
            female_num = ind_num_before-len(used_males)
        # filter female individuals according to missing rate or K-means cluster if no sample information
        else:
            if cutoff is None:
                female_num = 0
            else:
                del_list = []
                site_num = data.index.size
                missing_count = (lambda x: 'missing' if '.' in x or './.' in x else x)
                try:
                    # count misisng rate of each individual
                    if 'vcf' in file_format:
                        missing_rate_data = [data[i].map(missing_count).tolist().count('missing') / site_num for i in inds]
                    else:
                        missing_rate_data = [data[i].tolist().count('U') / site_num for i in inds]
                    # if set cutoff of missing rate, filter according to misisng rate
                    if cutoff:
                        del_list = [j for i, j in zip(missing_rate_data, inds) if i > cutoff]
                        female_num = len([i for i in inds if i not in del_list])
                    # if not set cutoff, will use K-means to cluster missing rate to identify samples and females
                    else:
                        del_list, female_num = self._get_female_by_KMeans(missing_rate_data, inds)
                    # if all individual filtered
                    if len(del_list) == ind_num_before:
                        self.logger.error('[Y-LineageTracker] Program stopped since no male sample left for subsequent analysis')
                        sys.exit()
                    # remove filtered indivisuals from data
                    data.drop(del_list, inplace=True, axis=1)
                except TypeError:
                    self.logger.error('[Y-LineageTracker] [Error] VCF file format error, please check your file')
                    sys.exit()

        inds_after = data.columns.tolist()[header_cut:]
        left_num = len(inds_after)
        if not cutoff is None:
            self.logger.info('[Y-LineageTracker] %d female individuals filtered and %d individuals left' % (female_num, left_num))

        return data, inds_after

    # main function for reading and filtering variant calling
    def read_filter_variant(self, cutoff, file_format):

        data_read, header_cut = self._read_variant(file_format)
        data_info, inds = self._filter_variant(data_read, cutoff, file_format, header_cut)

        return data_info, inds

    # filter sites and keep only SNPs for time estimation
    def restrict_time_region(self, data):

        # read data of region for time calculation
        from FilesIO import CommonData
        common_data = CommonData()
        region_info = common_data.read_region_info()
        time_data_info = pd.DataFrame(columns=data.columns.tolist())

        # function used to filter indels
        def remove_indels(x):

            alt = x.split(',')
            flag_num = sum([0 if len(i) == 1 else 1 for i in alt])
            if flag_num == 0:
                return True
            else:
                return False

        # restrict sites in regions
        for i in region_info.index:
            start_pos = int(region_info.at[i, 'start'])
            end_pos = int(region_info.at[i, 'end'])
            time_region_data = data[pd.to_numeric(data['POS']).map(lambda x: x >= start_pos and x <= end_pos)]
            time_data_info = pd.concat([time_data_info, time_region_data])

        # remove indels
        time_data_info = time_data_info[time_data_info['REF'].map(lambda x: len(x) == 1)]
        time_data_info = time_data_info[time_data_info['ALT'].map(remove_indels)]

        return time_data_info
