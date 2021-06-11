import sys

version = '1.3.0'
help_info = '\n' \
            '[Y-LineageTracker] a framework to fully analyze human Y-chromosome sequencing data from paternal level\n' \
            '\n' \
            'Version: ' + version + '\n' \
            'Usage: LineageTracker <subcommand> [options]\n' \
            '\n' \
            '[Subcommands]\n' \
            '--NRY haplogroup analysis\n' \
            '    classify    Classify NRY haplogroups and construct lineages from BAM or VCF.\n' \
            '    cluster     Perform NRY haplogroup clustering analysis.\n' \
            '    phylo       Perform Phylogenetic analysis for NRY haplogroups.\n' \
            '\n' \
            '--Y-STR and Y-haplotype analysis\n' \
            '    genostr     Genotype Y-STR from BAM or VCF indels.\n' \
            '    net         Perform network analysis from Y haplotype data.\n' \
            '    stat        Perform statistical analysis from Y haplotype data.\n' \
            '\n' \
            '--Time estimation\n' \
            '    time        Estimate divergence time of NRY haplogroups.\n' \
            '    tmrca       Estimate TMRCA of Y haplotypes.\n' \
            '\n' \
            'For more details about subcommands, type: LineageTracker <subcommand> -h/--help\n' \
            '\n' \
            '[General Options]\n' \
            '    -h/--help       Print help massages.\n' \
            '    -v/--version    Display the version of Y-LineageTracker you are using.\n'


def main():

    if len(sys.argv) > 1:
        # command for testing
        if sys.argv[1] == 'test':
            from Test import TestProgram
            TestProgram.main()
        # commands for analysis
        elif sys.argv[1] == 'classify':
            import ClassifyHaplogroup
            ClassifyHaplogroup.main()
        elif sys.argv[1] == 'cluster':
            import ClusterHaplogroup
            ClusterHaplogroup.main()
        elif sys.argv[1] == 'phylo':
            import PhyloHaplogroup
            PhyloHaplogroup.main()
        elif sys.argv[1] == 'genostr':
            import GenotypeSTR
            GenotypeSTR.main()
        elif sys.argv[1] == 'net':
            import HaploNetwork
            HaploNetwork.main()
        elif sys.argv[1] == 'stat':
            import StatisticsAnalysis
            StatisticsAnalysis.main()
        elif sys.argv[1] == 'time':
            import EstimateTime
            EstimateTime.main()
        elif sys.argv[1] == 'tmrca':
            import CommonAncestor
            CommonAncestor.main()
        elif sys.argv[1] == 'test':
            from Test import TestProgram
        # command for version information
        elif sys.argv[1] in ['version', '--version', '-v']:
            print('[Y-LineageTracker] version: %s' % version)
        # command for help information
        elif sys.argv[1] in ['help', '--help', '-h']:
            print(help_info)
        else:
            print(help_info)
    else:
        print(help_info)

if __name__ == '__main__':
    main()
