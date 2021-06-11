import os
import sys
import gzip
import tempfile
import numpy as np
import pandas as pd

sys.path.append(os.path.split(os.path.realpath(__file__))[0]+'/..')
from Bio import AlignIO, SeqIO, Phylo
from GetConfig import getConfig


config = getConfig()


class ConvertData(object):
    '''
    This class is mainly used to convert different format of data
    Conversion of data is helpful for some analysis
    This class can:
    1. convert VCF to MSA and matrix
    2. convert matrix to MSA
    3. convert MSA to matrix
    4. convert MSA in different formats
    MSA includes: fasta, phylip, nexus, mega
    VCF: variant call format
    MSA: multiple sequence alignment
    '''

    def __init__(self, in_file):

        from FilesIO import  check_compress
        self.vcf_header_cut = config.getint('DataInfo', 'VcfHeaderCut')
        self.compressed = check_compress(in_file)
        self.in_file = in_file

    # convert VCF to sequence format
    def _match_to_vcf(self, file, out_format):

        # dictionary of IUPAC ambiguities for nucleotides
        # '*' means deletion for GATK
        # deletions are ignored when making the consensus
        amb = {('A','A'):'A',
               ('A','N'):'A',
               ('C','C'):'C',
               ('C','N'):'C',
               ('G','G'):'G',
               ('G','N'):'G',
               ('N','A'):'A',
               ('N','C'):'C',
               ('N','G'):'G',
               ('N','N'):'N',
               ('N','T'):'T',
               ('T','N'):'T',
               ('T','T'):'T',
               ('*','*'):'-',
               ('A','*'):'A',
               ('*','A'):'A',
               ('C','*'):'C',
               ('*','C'):'C',
               ('G','*'):'G',
               ('*','G'):'G',
               ('T','*'):'T',
               ('*','T'):'T',
               ('N','*'):'N',
               ('*','N'):'N'}

        if self.compressed:
            vcf = gzip.open(file, 'rt')
        else:
            vcf = open(file, 'r')
        with vcf:

            # create a list to store sample names
            sample_names = []

            # keep track of longest sequence name for padding with spaces in the output file
            len_longest_name = 0

            # look for the line in the VCF header with the sample names
            for line in vcf:
                if line.startswith('#CHROM'):

                    # split line into fields
                    broken = line.strip('\n').split('\t')
                    min_samples_locus = int(len(broken[self.vcf_header_cut:])*0.8)

                # Create a list of sample names and the keep track of the longest name length
                    for i in range(self.vcf_header_cut, len(broken)):
                        name_sample = broken[i].replace('./', '') # GATK adds './' to sample names sometimes
                        sample_names.append(name_sample)
                        len_longest_name = max(len_longest_name, len(name_sample))

                # Find out the ploidy of the genotypes
                elif not line.startswith('#'):
                    break
        vcf.close()

        # create an intermediate file to hold the sequence data vertically and then transpose it to create the matrices
        if out_format == 'matrix':
            tmp_list = []
        else:
            tmp_file = tempfile.NamedTemporaryFile(mode='w')
        # process vcf
        index_last_sample = len(sample_names) + self.vcf_header_cut
        # initialize line counter
        snp_num = 0
        snp_accepted = 0
        snp_shallow = 0
        snp_multinuc = 0
        snp_biallelic = 0
        ploidy = 1
        if self.compressed:
            vcf = gzip.open(file, 'rt')
        else:
            vcf = open(file, 'r')
        with vcf:
            while True:
                # load large chunks of file into memory
                vcf_chunk = vcf.readlines(100000)
                if not vcf_chunk:
                    break
                # process the SNPs one by one
                for line in vcf_chunk:
                    if not line.startswith('#') and line.strip('\n') != '': # pyrad sometimes produces an empty line after the #CHROM line
                        # split line into columns
                        broken = line.strip('\n').split('\t')
                        for g in range(self.vcf_header_cut, len(broken)):
                            if broken[g] in ['.', '.|.']:
                                broken[g] = './.'
                        # keep track of number of genotypes processed
                        snp_num += 1
                        # check if the SNP has the minimum of samples required
                        if (len(broken[self.vcf_header_cut:]) - ''.join(broken[self.vcf_header_cut:]).count('./.')) >= min_samples_locus:
                            # check that ref genotype is a single nucleotide and alternative genotypes are single nucleotides
                            if len(broken[3]) == 1 and (len(broken[4])-broken[4].count(',')) == (broken[4].count(',')+1):
                                # add to running sum of accepted SNPs
                                snp_accepted += 1
                                # create a dictionary for genotype to nucleotide translation
                                # each SNP may code the nucleotides in a different manner
                                nuc = {str(0): broken[3], '.': 'N'}
                                for n in range(len(broken[4].split(','))):
                                    nuc[str(n+1)] = broken[4].split(',')[n]
                                # translate genotypes into nucleotides and the obtain the IUPAC ambiguity
                                # for heterozygous SNPs, and append to DNA sequence of each sample
                                if ploidy == 1:
                                    site_list = [(amb[(nuc[broken[i][0]], nuc[broken[i][0]])]) for i in range(self.vcf_header_cut, index_last_sample)]
                                    site_tmp = ''.join(site_list)
                                # write entire row of single nucleotide genotypes to temporary file of fasta and nexus
                                if out_format == 'matrix':
                                    test_list.append(site_list)
                                else:
                                    tmp_file.write(site_tmp+'\n')
                            else:
                                # keep track of loci rejected due to multinucleotide genotypes
                                snp_multinuc += 1
                                # keep track of loci rejected due to exceeded missing data
                                snp_shallow += 1
                        else:
                            # keep track of loci rejected due to exceeded missing data
                            snp_shallow += 1

        if out_format == 'matrix':
            matrix = np.array(test_list).T
            matrix_df = pd.DataFrame(matrix, index=sample_names, columns=['V'+str(i) for i in range(np.size(test2, 1))])
            matrix_df.index.name = 'SampleID'
        else:
            return tmp_file, sample_names, len_longest_name, snp_accepted

    # convert VCF to sequence format
    def _write_matched_tmp_to_seq(self, tmp_file, sample_names, len_longest_name, snp_accepted, out_format):

        tmp_file2 = tempfile.NamedTemporaryFile(mode='w')

        # output file type for sebsequent analysis
        if out_format == 'phylip':
            header_phy = str(len(sample_names))+' '+str(snp_accepted)+'\n'
            tmp_file2.write(header_phy)
        elif out_format == 'nexus':
            header_nex = '#NEXUS\n\nBEGIN DATA;\n\tDIMENSIONS NTAX=' + str(len(sample_names)) + ' NCHAR=' + str(snp_accepted) + ';\n\tFORMAT DATATYPE=DNA' + ' MISSING=N' + ' GAP=- ;\nMATRIX\n'
            tmp_file2.write(header_nex)

        # Write sequences of the ingroup
        for s in range(0, len(sample_names)):
            with open(tmp_file.name) as tmp_seq:
                seqout = ''
                # where the transposing happens
                for line in tmp_seq:
                    seqout += line[s]
                if out_format == 'fasta':
                    # write FASTA line
                    tmp_file2.write('>'+sample_names[s]+'\n'+seqout+'\n')
                else:
                    # pad sequences names and write PHYLIP or NEXUS lines
                    padding = (len_longest_name+3-len(sample_names[s])) * ' '
                    if out_format == 'phylip':
                        tmp_file2.write(sample_names[s]+padding+seqout+'\n')
                    if out_format == 'nexus':
                        tmp_file2.write(sample_names[s]+padding+seqout+'\n')

        res = AlignIO.read(open(tmp_file2.name), out_format)

        return res

    # convert VCF to msa
    def _vcf_to_msa(self, file, out_format):

        tmp_file, sample_names, len_longest_name, snp_accepted = self._match_to_vcf(file, out_format)
        res = self._write_matched_tmp_to_seq(tmp_file, sample_names, len_longest_name, snp_accepted, out_format)

        return res

    # convert VCF to matrix format
    def _vcf_to_matrix(self, file, out_format='matrix'):

        matrix_df = self._match_to_vcf(file, out_format)

        return matrix_df

    # convert msa to matrix
    def _seq_to_matrix(self, file, in_format, tmp):

        if tmp:
            msa = AlignIO.read(open(file.name), in_format)
        else:
            msa = AlignIO.read(open(file), in_format)
        ids = [i.id for i in msa]
        seqs = [i.seq for i in msa]
        seq_len = len(seqs[0])
        seq_matrix = np.array([[j]+[x for x in i] for i, j in zip(seqs, ids)])

        matirx = pd.DataFrame(columns=['SampleID']+['V'+str(i) for i in range(seq_len)], data=seq_matrix)

        return matirx

    # read seq to msa
    def _seq_to_msa(self, file, in_format, tmp):

        if tmp:
            res = AlignIO.read(open(file.name), in_format)
        else:
            res = AlignIO.read(open(file), in_format)

        return res

    # convert msa in different
    def _seq_convert(self, file, in_format, out_format, tmp):

        tmp_file = tempfile.NamedTemporaryFile(mode='w')
        if tmp:
            SeqIO.convert(file.name, in_format+'-relaxed', tmp_file, out_format)
        else:
            SeqIO.convert(file, in_format+'-relaxed', tmp_file, out_format)
        res = AlignIO.read(open(tmp_file.name), out_format)

        return res

    # convert matrix to msa
    def _matrix_to_msa(self, file, out_format):

        matrix_info = np.loadtxt(file, dtype='object').T
        snp_accepted = len(matrix_info[:,0]) - 1
        sample_names= matrix_info[0,:].tolist()[1:]
        tmp_file = tempfile.NamedTemporaryFile(mode='w')
        for i in matrix_info[1:, 1:]:
            tmp_file.write(''.join(i)+'\n')

        res = self._write_matched_tmp_to_seq(tmp_file, sample_names, snp_accepted, out_format)

        return res

    # convert meg to fasta
    def _meg_to_fasta(self, file):

        meg = open(file, 'r').readlines()
        tmp_file = tempfile.NamedTemporaryFile(mode='w')

        with open(tmp_file.name, 'w') as tmp:
            for i in meg[3:]:
                if i.startswith('#'):
                    tmp.write(i.replace('#', '>', 1))
                else:
                    tmp.write(i)

        return tmp_file

    # main function of data conversion
    def data_convert(self, in_format, out_format):

        seq_format = ['fasta', 'phylip', 'nexus']
        in_file = self.in_file

        if in_format == 'vcf':
            if out_format == 'msa':
                res = self._vcf_to_msa(in_file, 'fasta')
            elif out_format == 'matrix':
                res = self._vcf_to_martix(in_file)
        elif in_format == 'matrix':
            res = self._matrix_to_msa(in_file, out_format)
        else:
            if in_format == 'meg':
                in_file = self._meg_to_fasta(in_file)
                in_format = 'fasta'
                tmp = True
            else:
                tmp = False
            if out_format == 'matrix':
                res = self._seq_to_matrix(in_file, in_format, tmp)
            elif out_format == 'msa':
                res = self._seq_to_msa(in_file, in_format, tmp)
            else:
                res = self._seq_convert(in_file, in_format, out_format, tmp)

        return res
