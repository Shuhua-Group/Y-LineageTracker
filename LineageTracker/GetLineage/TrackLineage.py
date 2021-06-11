import re
import sys
import logging
import numpy as np
import pandas as pd

from functools import reduce
from itertools import filterfalse
from collections import Counter, OrderedDict
from GetConfig import getConfig


config = getConfig()


class TrackLineage(object):
    '''
    This class is used to:
    1. classify NRY haplogroup
    2. Track lineage
    3. Output result of classification
    '''

    def __init__(self, matched_info, ref_ind_info, approximate_haps, build, format, chip):

        # record running log
        self.logger = logging.getLogger()

        # get reference haplogroups of which all indivisuals are same
        if build == '37' or 'Ref':
            idx = -7
        elif build == '38':
            idx = -8

        # get parameters from config file
        self.build = build
        self.marker_rate = float(config.get('HaplogroupTree', 'MarkerRate'))
        self.length_cutoff = config.getint('HaplogroupTree', 'LengthCutoff')
        self.matched_marker = config.get('HaplogroupMarker', 'MatchMarker')
        self.match_miss_marker = config.get('HaplogroupMarker', 'MatchMissMarker')
        self.track_symbol = config.get('HaplogroupMarker', 'TrackSymbol') # track symbol is '->'
        self.track_rate = float(config.get('HaplogroupTree', 'TrackRate'))
        self.tree_trunk = config.get('HaplogroupTree', 'CommonTrunk').split(',')
        self.approximate_haps = approximate_haps
        self.haplogroup_interval = config.getint('HaplogroupTree', 'HaplogroupInterval')

        if format == 'vcf' and not chip and build != 'Ref':
            common_sites = [i for i in ref_ind_info['Build'+str(build)].tolist() if i not in matched_info[:, idx].tolist()]
            common_df = ref_ind_info[ref_ind_info['Build'+str(build)].map(lambda x: x in common_sites)]
            self.common_haps = list(set(common_df['Haplogroup'])&set(self.tree_trunk))
            self.common_key_haps = list(set(common_df[common_df['KeyInfo']=='Key']['Haplogroup'])&set(self.tree_trunk))
        else:
            self.common_key_haps = []
            self.common_haps = []

    # extract matched haplogroups and their frequence, individuals with no matched item will be removed and recoreded
    def _extract_haplogroup(self, col, key_matched_info, ref_haps):

        # extract all matched key haplogroups
        key_haps = list(map(lambda x: x[1:], list(filter(lambda x: x[0] == self.matched_marker, key_matched_info[:, col])))) + self.common_key_haps
        count_res = OrderedDict(sorted(Counter(all_matched_haps+self.common_haps).items()))
        matched_haps = list(count_res.keys())
        matched_haps_freq = list(count_res.values())
        total_haps_freq = list(dict(filter(lambda x: x[0] in set(matched_haps), ref_haps)).values())

        unzipped_haps = [matched_haps, matched_haps_freq, total_haps_freq]

        return key_haps, unzipped_haps

    # find the most downstream haplogroup according to key haplogroups
    def _search_final_hap(self, key_haps):

        # tracker rate of haplogroup A is special since the complex names
        haplogroup_A_rate = float(config.get('HaplogroupTree', 'HaplogroupARate'))

        # set initial haplogroup is NONE
        key_hap = self.match_miss_marker
        final_hap = ''

        # traverse to find a most reliable key haplogroup
        length = 0
        for hap in filterfalse(lambda x: x.startswith('K') or x in self.tree_trunk, key_haps):
            length2 = len(hap)
            if length2 == 1:
                continue
            if length2 >= length:
                exist_num = 0
                hap_num = length2 - 1
                for i in hap:
                    hap_str = hap[0: hap_num]
                    if hap_str in key_haps:
                        exist_num += 1
                    hap_num -= 1
                if exist_num == length2-1:
                    length = length2
                    final_hap = hap
                    continue
                if length2 >= self.length_cutoff and (exist_num+1)/length2 > self.track_rate:
                    length = length2
                    final_hap = hap
                    continue
                if hap.startswith('A'):
                    if (exist_num+1)/length2 > (self.track_rate - haplogroup_A_rate):
                        length = length2
                        final_hap = hap
                    elif hap == 'A0' and 'A0T' in key_haps:
                        final_hap = hap
                        length = length2
                    elif hap.startswith('A1b') and 'A1b' in key_haps and 'A1' in key_haps:
                        final_hap = hap
                        length = length2

        return key_hap, final_hap

    # search whether there are more downstream haplogroup
    def _search_downstream(self, matched_haps, zipped_haps, key_haps, key_hap, final_hap):

        key_hap = final_hap
        length = len(final_hap)
        rate = 0
        downstream = sorted(filter(lambda x: len(x[0]) >= len(final_hap)+1 and
                                                 x[0].startswith(final_hap) and
                                                 x[0] not in self.tree_trunk, zipped_haps),
                                                 key=lambda x: len(x[0]))

        traversed = []
        for hap, matched, total in downstream:
            length = len(hap)
            final_hap = hap
            rate = matched/total
            if final_hap in key_haps and len(hap)-length < self.haplogroup_interval:
                key_hap = final_hap
            else:
                if matched/total > rate and len(hap)-length < self.haplogroup_interval:
                    length = len(hap)
                    final_hap = hap
                    rate = matched/total
                    if final_hap in key_haps:
                        key_hap = final_hap
                elif len(list(filter(lambda x: hap.startswith(x), traversed))) > len(list(filter(lambda x: final_hap.startswith(x), traversed)))-1: # exclude final hap itself
                    length = len(hap)
                    final_hap = hap
                    rate = matched/total
                    if final_hap in key_haps:
                        key_hap = final_hap
            traversed.append(hap)

        return key_hap, final_hap

    # if the final haplogroup is in K branch
    def _replace_single_K_trunk_haplogroup(self, key_haps, key_hap):

        if key_hap == 'K2b1' or key_hap == 'K2b':
            if 'S' in key_haps:
                key_hap = 'S'
            elif 'M' in key_haps:
                key_hap = 'M'
            elif 'P1' in key_haps:
                key_hap = 'P1'
            elif 'P' in key_haps:
                key_hap = 'P'
        elif key_hap == 'K2':
            if 'NO1' in key_haps:
                key_hap = 'NO1'
            elif 'NO' in key_haps:
                key_hap = 'NO'

        final_hap = key_hap

        return key_hap, final_hap

    # if the resolution of data is low, detect potential haplogroups
    def _detect_low_resolution_haplogroup(self, zipped_haps, matched_haps, key_haps, final_hap, key_hap='.'):

        length = 0
        rate = 0
        flag = False
        downstream = sorted(zipped_haps, key=lambda zipped_haps: len(zipped_haps[0]))
        downstream = filter(lambda x: x[0] not in self.tree_trunk, downstream)
        for hap, matched, total in downstream:
            if hap[0] not in matched_haps:
                continue
            if matched/total <= (self.marker_rate-0.45):
                continue
            if hap.startswith('A'):
                matched_haps.append('A')
            length2 = len(hap)
            if (length2 > length) or (length2 == length and matched/total > rate):
                exist_num = 0
                hap_num = length2 - 1
                for i in hap:
                    hap_str = hap[0: hap_num]
                    if hap_str in matched_haps:
                        exist_num += 1
                    hap_num -= 1
                if exist_num == length2-1:
                    flag = True
                else:
                    if length2 >= 4 and (exist_num+1)/length2 > (self.track_rate-0.2):
                        flag  = True
                if flag:
                    length = length2
                    final_hap = hap
                    rate = matched/total
                    flag = False
                    for j in reversed(sorted(key_haps)):
                        if final_hap.startswith(j):
                            key_hap = j
                            break

        if final_hap:
            key_hap, final_hap = self._search_downstream(matched_haps, zipped_haps, key_haps, key_hap, final_hap)

        return key_hap, final_hap

    # if the resolution of data is lower than above, detect potential haplogroups
    def _detect_lower_resolution_haplogroup(self, zipped_haps, key_haps, key_hap):

        if len(list(filterfalse(lambda x: x in self.tree_trunk, key_haps))) > 1:
            length = 0
            last_haps = list(filterfalse(lambda x: x in self.tree_trunk, key_haps))
            for hap in last_haps:
                length2 = len(hap)
                if length2 >= length:
                    final_hap = hap
                    length = length2
                    key_hap = self.matched_marker + final_hap
        else:
            hap_rate = 0
            last_haps = sorted(zipped_haps, key=lambda zipped_haps: zipped_haps[0])
            last_haps = sorted(last_haps, key=lambda last_haps: len(zipped_haps[0]))
            if set(list(zip(*last_haps))[0]).issubset(set(self.tree_trunk)):
                final_hap = last_haps[-1][0]
            else:
                for hap, matched, total in last_haps:
                    if len(hap) > 1 and hap in self.tree_trunk:
                        continue
                    if matched/total >= hap_rate:
                        final_hap = hap
                        hap_rate = matched/total

        return key_hap, final_hap

    # get a lineage track according to ISOGG tree information
    def _track_lineage(self, final_hap):

        from FilesIO import CommonData
        common_data = CommonData()
        tree_info = common_data.read_tree_info()

        if (final_hap.isalpha() and final_hap.isupper()) or final_hap == 'NO1':
            track = final_hap
        else:
            track = ''
            alpha = not final_hap[-1].isalpha()
            for num_str in range(len(final_hap)):
                tmp_up = final_hap[0: len(final_hap)-num_str]
                is_alpha = tmp_up[-1].isalpha()
                if is_alpha != alpha:
                    up = tmp_up
                else:
                    continue
                if up in self.approximate_haps:
                    up = up + '~'
                track = track + self.track_symbol + up
                alpha = is_alpha
            track = track[2:]
            # A is special haplogroup in name
            if track[-1] == config.get('HaplogroupA', 'SpecificHaplogroup'):
                track_temp = track[:-3]
                track = self.track_symbol + track_temp + config.get('HaplogroupA', 'SpecificTrackOne')
                pattern_track = re.compile(r'A00[abc]')
                if re.search(pattern_track, track):
                    track_temp = pattern_track.search(track).group()
                    track = self.track_symbol + track_temp + config.get('HaplogroupA', 'SpecificTrackTwo')
                if config.get('HaplogroupA', 'SpecificSubhaplogroup') in track:
                    track = config.get('HaplogroupA', 'SpecificTrackThree')

        pattern = re.compile(r'[A-Za-z0-9]+')
        length = 0
        num = 0
        header = track.split('>')[-1]
        for i in tree_info:
            if len(i) < len(header)+1:
                continue
            if i[-(len(header)+1)] == '\t' and header == i[-len(header):]:
                    length = i.count('-')
                    for j in list(reversed(tree_info)):
                        if j.count('-') + 1 == length and tree_info.index(j) < tree_info.index(i):
                            searched = pattern.search(j).group()
                            track = track + self.track_symbol + searched
                            length = length - 1
                        else:
                            continue
                    track = self.track_symbol + track + self.track_symbol + 'Y-Adam'

        return track

    # count how many mutations matched
    def _count_matched_number(self, track, unzipped_haps):

        matched_haps, matched_haps_freq, total_haps_freq = unzipped_haps
        haps_in_track = track.split(self.track_symbol)[:-1]
        for hap in haps_in_track:
            if hap in matched_haps:
                hap_index = matched_haps.index(hap)
            elif hap.strip('~') in matched_haps:
                hap_index = matched_haps.index(hap.strip('~'))
            else:
                continue
            matched = matched_haps_freq[hap_index]
            total = total_haps_freq[hap_index]
            hap_without_count = self.track_symbol + hap + self.track_symbol
            hap_with_count = self.track_symbol + hap + '(' + str(matched) + '/' + str(total) + ')' + self.track_symbol
            track = track.replace(hap_without_count, hap_with_count)
        track = track[2:]

        return track

    # main function of classification: classify haplogroup of each individual
    def classify_haplogroup(self, inds, matched_info, ref_info, output_set, mutation, header_cut):

        # output data one: haplogroup classification result
        hap_data = pd.DataFrame(index=inds,
                                columns=['Haplogroup'])
        output_hap_file= output_set[0]

        # output data two: lineage of haplogroups
        lineage_data = pd.DataFrame(index=inds,
                                    columns=['Haplogroup',
                                            'KeyHaplogroup',
                                            'Mutations',
                                            'LineageTrack'])
        output_lineage_file = output_set[1]
        # output data three: all matched mutation
        if mutation:
            mutation_data = pd.DataFrame(index=inds,
                                         columns=['AllMutations'])
            output_mutation_file = output_set[2]

        # initial variables
        no_hap_inds = []

        # get freq of all total matched haplogroups
        ref_haps = sorted(Counter(ref_info['Haplogroup']).items())

        base_cols = list(range(0, 9)) + list(range(-11, 0))
        base_info = matched_info[:, base_cols]
        key_matched_info = matched_info[matched_info[:,-4]=='Key']

        # traverse all indivisual to classify haplogroups
        final_list = []
        key_list = []
        mutation_list = []
        track_list = []
        all_mutation_list = []
        ind_num = 0
        for n, ind in enumerate(inds):

            col = n + header_cut
            col_info = matched_info[:, col]

            # extract all matched haplogroups
            global all_matched_haps
            # [1:] means remove match marker
            all_matched_haps = sorted(list(map(lambda x: x[1:], list(filter(lambda x: x[0]==self.matched_marker, col_info)))))

            # whether there are individual who don't have matched haplogroups
            if len(all_matched_haps) == 0:
                no_hap_inds.append(ind)
                final_list.append('.')
                key_list.append('.')
                mutation_list.append('.')
                track_list.append('.')
                self.logger.warning('[Y-LineageTracker] [Warning] There are no mutation matched for sample %s' % (ind))
                continue

            # extract matched haplogroups and their frequence
            key_haps, unzipped_haps = self._extract_haplogroup(col, key_matched_info, ref_haps)
            matched_haps, matched_haps_freq, total_haps_freq = unzipped_haps

            zipped_haps = list(zip(*unzipped_haps))

            if len(key_haps) > 0: # if there are key haplogroups
                key_hap, final_hap = self._search_final_hap(key_haps) # whether there is a key haplogroup, data with low resolution will get an empty string
                if final_hap: # haplogroup with single character will keep an empty string
                    key_hap, final_hap = self._search_downstream(matched_haps, zipped_haps, key_haps, key_hap, final_hap) # whether the final haplogroup have a more donwstream one
                else: # there are key haplogroups but not classified
                    key_hap, final_hap = self._detect_low_resolution_haplogroup(zipped_haps, matched_haps, key_haps, final_hap, key_hap)
            else: # no key haplogroup matched
                if self.build != 'Ref':
                    self.logger.warning('[Y-LineageTracker] [Warning] There are no key haplogroups matched for sample %s' % (ind))
                final_hap = ''
                key_hap, final_hap = self._detect_low_resolution_haplogroup(zipped_haps, matched_haps, key_haps, final_hap)

            if key_hap.startswith('K'): # other haplogroups not started with K have been tested
                key_hap, final_hap = self._replace_single_K_trunk_haplogroup(key_haps, key_hap)

            if final_hap == '': # if the haplogroup is None
                key_hap, final_hap = self._detect_lower_resolution_haplogroup(zipped_haps, key_haps, key_hap)

            # match all matched haplogroups to their mutations
            matched_haps_array = base_info[list(map(lambda x: x[0] == self.matched_marker, col_info))]
            mutations = ', '.join(matched_haps_array[:, -11])

            # match final haplogroup to its mutations
            final_hap_array = matched_haps_array[matched_haps_array[:, -10] == final_hap]
            final_mutations = ', '.join(final_hap_array[:, -11])
            if 'Key' not in final_hap_array[:, -4]:
                final_mutations = '*' + final_mutations

            # get lineage track of haplogroup
            track = self._track_lineage(final_hap)

            # count matched number of each haplogroup
            haps_in_track = track.split(self.track_symbol)[:-1]
            track = self._count_matched_number(track, unzipped_haps)

            # output row for each individual
            mutation_list.append(final_mutations)
            track_list.append(track)
            if mutation:
                all_mutation_list.append(mutations)
            if final_hap in self.approximate_haps:
                final_hap = final_hap + '~'

            final_list.append(final_hap)
            key_list.append(key_hap)

            percent = float(ind_num+1)*100 / float(len(inds))
            sys.stdout.write('[Y-LineageTracker] Classifying haplogroups…… %.2f' % percent)
            sys.stdout.write('%\r')
            sys.stdout.flush()
            ind_num += 1

        hap_data['Haplogroup'] = final_list
        # output file one: haplogroup classification result
        hap_data.index.name = 'SampleID' # set individuals as index
        hap_data.to_csv(output_hap_file, sep='\t', index=True, encoding='utf-8')

        print()
        self.logger.info('[Y-LineageTracker] Classification and Tracking haplogroup lineages finished')

        # output file two: lineage of each haplogroup
        lineage_data['Haplogroup'] = final_list
        lineage_data['KeyHaplogroup'] = key_list
        lineage_data['Mutations'] = mutation_list
        lineage_data['LineageTrack'] = track_list
        lineage_data.index.name = 'SampleID'
        lineage_data.to_csv(output_lineage_file, sep='\t', index=True)
        # output file three: matched mutations
        if mutation:
            mutation_data.index.name = 'SampleID'
            mutation_data['AllMutations'] = all_mutation_list
            mutation_data.to_csv(output_mutation_file, sep='\t', index=True)

        if len(no_hap_inds) != 0:
            if len(no_hap_inds) <= 5:
                self.logger.info('[Y-LineageTracker] [Warning] %s cannot classify haplogroup due to few or mismatched sites' % ', '.join(no_hap_inds))
            else:
                self.logger.info('[Y-LineageTracker] [Warning] %d samples cannot classify haplogroup due to few or mismatched sites' % len(no_hap_inds))

        self.logger.info('[Y-LineageTracker] Haplogroup classification done')

        hap_data.reset_index(level=0, inplace=True)

        return hap_data
