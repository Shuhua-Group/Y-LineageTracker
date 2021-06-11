import re
import sys
import time
import pysam
import logging
import pandas as pd
import numpy as np

from GetConfig import getConfig


config = getConfig()


class MatchVariant(object):
    '''
    This class is used to:
    1. for NRY haplogroup: match sites and variants of input data to ISOGG reference
    2. for Y-STRs: match VCF indels of input data to Y-STR panel
    '''

    def __init__(self, build):

        self.logger = logging.getLogger()
        self.build = build

    # match ISOGG reference to vcf, and count haplogroup of each individual
    def match_variant_snp(self, ref_info, data_info, file_format):

        # match common position of ISOGG haplogroup data and input data
        if file_format == 'bam':
            header_cut = config.getint('DataInfo', 'BamHeaderCut')
            matched_info = data_info
        else:
            matched_marker = config.get('HaplogroupMarker', 'MatchMarker') # used to mark matched haplogroups

            matched_info = pd.merge(data_info, ref_info,
                                    left_on='POS',
                                    right_on='Build'+str(self.build),
                                    how='inner',
                                    sort=False)
            matched_info = np.array(matched_info)
            header_end = -config.getint('DataInfo', 'ISOGGInfoCut')
            header_cut = config.getint('DataInfo', 'VcfHeaderCut')
            hap_function = (lambda x: matched_marker+matched_info[i, -10] if x[0] == match_marker else x)
            ins_function = (lambda x: matched_marker+matched_info[i, -10] if x[0] in ins_loc_list else x)
            del_function = (lambda x: matched_marker+matched_info[i, -10] if x[0] in del_loc_list else x)

            for i in range(len(matched_info)):
                if ('ins' or 'del') not in matched_info[i, -6]:
                    if matched_info[i, -6] == matched_info[i, 3]:
                        match_marker = '0'
                        matched_info[i, header_cut: header_end] = list(map(hap_function, matched_info[i, header_cut: header_end]))
                    else:
                        variant_list = matched_info[i, 4].split(',')
                        if len(variant_list) > 1:
                            if matched_info[i, -6] in variant_list:
                                alt_num = str(variant_list.index(matched_info[i, -6]) + 1)
                                match_marker = alt_num
                                matched_info[i, header_cut: header_end] = list(map(hap_function, matched_info[i, header_cut: header_end]))
                        elif len(variant_list) == 1:
                            if matched_info[i, -6] == matched_info[i, 4]:
                                match_marker = '1'
                                matched_info[i, header_cut: header_end] = list(map(hap_function, matched_info[i, header_cut: header_end]))

                elif matched_info[i, -6].startswith('ins'):
                    if matched_info[i, -6] == 'ins':
                        variant_list = matched_info[i, 4].split(',')
                        if len(variant_list) > 1:
                            ins_list = list(map(lambda x: len(matched_info[i, 3]) < len(x), variant_list))
                            ins_loc_list = np.argwhere(ins_list == True).shape
                            if True in ins_list:
                                matched_info[i, header_cut: header_end] = list(map(ins_function, matched_info[i, header_cut: header_end]))
                        elif len(variant_list) == 1:
                            if len(matched_info[i, 4]) > len(matched_info[i, 3]):
                                match_marker = '1'
                                matched_info[i, header_cut: header_end] = list(map(hap_function, matched_info[i, header_cut: header_end]))
                    elif re.search(r'ins \d+ bp', matched_info[i, -6]):
                        ins_bp = matched_info[i, -6][4]
                        variant_list = matched_info[i, 4].split(',')
                        if len(variant_list) > 1:
                            ins_list = list(map(lambda x: str(len(x) - len(matched_info[i, 3])), variant_list))
                            ins_loc_list = np.argwhere(ins_list == ins_bp).shape
                            if ins_bp in ins_list:
                                matched_info[i, header_cut: header_end] = list(map(ins_function, matched_info[i, header_cut: header_end]))
                        elif len(variant_list) == 1:
                            if len(matched_info[i, 4]) - len(matched_info[i, 3]) == ins_bp:
                                match_marker = '1'
                                matched_info[i, header_cut: header_end] = list(map(hap_function, matched_info[i, header_cut: header_end]))

                elif matched_info[i, -6].startswith('del'):
                    if matched_info[i, -6] == 'del':
                        variant_list = matched_info[i, 4].split(',')
                        if len(variant_list) > 1:
                            del_list = list(map(lambda x: len(matched_info[i, 3]) > len(x), variant_list))
                            del_loc_list = np.argwhere(del_list == True).shape
                            if True in del_list:
                                matched_info[i, header_cut: header_end] = list(map(del_function, matched_info[i, header_cut: header_end]))
                        elif len(variant_list) == 1:
                            if len(matched_info[i, 4]) < len(matched_info.at[i, 3]):
                                match_marker = '1'
                                matched_info[i, header_cut: header_end] = list(map(hap_function, matched_info[i, header_cut: header_end]))

                    elif re.search(r'del -\d+ bp', matched_info[i, -6]):
                            del_bp = matched_info[i, -6][5]
                            variant_list = matched_info[i, 4].split(',')
                            if len(variant_list) > 1:
                                del_list = list(map(lambda x: str(len(matched_info[i, 3]) - len(x)), variant_list))
                                del_loc_list = np.argwhere(del_list == del_bp).shape
                                if del_bp in del_list:
                                    matched_info[i, header_cut: header_end] = list(map(del_function, matched_info[i, header_cut: header_end]))
                            elif len(variant_list) == 1:
                                if len(matched_info[i, 3]) - len(matched_info[i, 4]) == del_bp:
                                    match_marker = '1'
                                    matched_info[i, header_cut: header_end] = list(map(hap_function, matched_info[i, header_cut: header_end]))

                percent = float(i + 1) * 100 / float(len(matched_info))
                sys.stdout.write('[Y-LineageTracker] Matching variants and sites to ISOGG reference…… %.2f' % percent)
                sys.stdout.write('%\r')
                sys.stdout.flush()

        print()
        self.logger.info('[Y-LineageTracker] Matching variants and sites to ISOGG reference finished')
        time.sleep(0.01)

        return matched_info, header_cut

    def match_variant_str(self, data_info, inds, panel):

        # read reference Y-STR information
        from FilesIO import get_str_from_panel
        str_info = get_str_from_panel(panel)

        to_merge_pos = []
        pre_start_pos = str_info['Start'+str(self.build)].tolist()
        motif_length = str_info['MotifLength'].tolist()
        start_pos = [i-j for i, j in zip(pre_start_pos, motif_length)]
        end_pos = str_info['End'+str(self.build)].tolist()
        [to_merge_pos.extend(list(range(i, j+1))) for i, j in zip(start_pos, end_pos)]
        to_merge_pos = [str(i) for i in to_merge_pos]
        pos_set = sorted(list(set(set(data_info['POS'].tolist()) & set(to_merge_pos))))

        data_info = data_info[data_info['POS'].map(lambda x: x in pos_set)]

        data_info['POS'] = pd.to_numeric(data_info['POS'])

        str_haplotype = pd.DataFrame(columns=inds)
        # get VCF header number and create a empty dataframe for output
        header_cut = config.getint('DataInfo', 'VcfHeaderCut')
        # traverse all the Y-STR
        for num, STR_idx in enumerate(str_info.index):
            # region restriction (to determine which should be used)
            start = str_info.at[STR_idx, 'Start'+str(self.build)]
            end = str_info.at[STR_idx, 'End'+str(self.build)]
            motif = str_info.at[STR_idx, 'RefMotif'].split(', ')
            motif_num = len(motif)
            str_name = str_info.at[STR_idx, 'Name'+str(self.build)]
            if motif_num == 1:
                motif = motif[0]
                minus_length = len(motif)
                STR = data_info[data_info['POS'].map(lambda x: x>=(start-minus_length) and x<=end)]
            else:
                minus_length = int(np.mean([len(i) for i in motif]))
                STR = data_info[data_info['POS'].map(lambda x: x>=(start-minus_length) and x<=end)]

            # match str allele
            base_motif_num = str_info.at[STR_idx, 'Ref'+str(self.build)]
            if STR.empty:
                motif_count = [base_motif_num for i in range(STR.columns.size-header_cut)]
                str_name = '*'+str_name
            else:
                if motif_num == 1:
                    filter_boole = [motif in i or motif in j for i, j in zip(STR['REF'], STR['ALT'])]
                else:
                    for single_motif_num in range(len(motif)):
                        single_motif = motif[single_motif_num]
                        single_boole = [single_motif in i or single_motif in j for i, j in zip(STR['REF'], STR['ALT'])]
                        if single_motif_num > 0:
                            filter_boole = [i or j for i, j in zip(single_boole, filter_boole)]
                        else:
                            filter_boole = single_boole
                    STR = STR[filter_boole]
                for STR_num, STR_index in enumerate(STR.index):
                    alt_list = STR.at[STR_index, 'ALT'].split(',')
                    if motif_num == 1:
                        ref_motif_num = STR.at[STR_index, 'REF'].count(motif)
                        alt_count = [i.count(motif)-ref_motif_num if i != '*' else 0 for i in alt_list]
                    else:
                        ref_motif_num_list  = [STR.at[STR_index, 'REF'].count(i) for i in motif]
                        alt_count_list = [[i.count(m)-n if i != '*' else 0 for i in alt_list] for m, n in zip(motif, ref_motif_num_list)]
                        alt_count = [sum(i) for i in list(zip(*alt_count_list))]
                    alt_count.insert(0, 0) # add REF as 0 alt
                    if STR_num == 0:
                        alt_num = ['.' if i.startswith('.') else alt_count[int(i[0])] for i in STR.iloc[STR_num, header_cut:]]
                    else:
                        alt_count_match = ['.' if i.startswith('.') else alt_count[int(i[0])] for i in STR.iloc[STR_num, header_cut:]]
                        alt_num = [sum([k for k in [i, j] if isinstance(k, int)]) if any (i.isdigit() for i in ''.join([str(i), str(j)])) else '.' for i, j in zip(alt_num, alt_count_match)]
                motif_count = [i+base_motif_num if not i == '.' else '.' for i in alt_num]
            str_haplotype.loc[str_name] = motif_count
            percent = float(num + 1) * 100 / float(str_info.index.size)
            sys.stdout.write('[Y-LineageTracker] Genotyping Y-STR…… %.2f' % percent)
            sys.stdout.write('%\r')
            sys.stdout.flush()


        str_haplotype = str_haplotype.T

        print()
        self.logger.info('[Y-LineageTracker] Genotyping Y-STR finished')

        return str_haplotype
