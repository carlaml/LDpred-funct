#!/usr/bin/env python
# Usar ldpred-master/coord_genotypes_ldpredfunct.py como base

#### bim_file = "/n/scratch2/sg374/UKB_LIKELIHOOD/25K_Carla/plink_files/Final_UK10K.[1:22]"

### input plink files: per chromsom
import getopt
import sys
import os
import h5py
import scipy as sp
from scipy import stats
from plinkio import plinkfile
import itertools as it
import gzip
import glob
import re

ambig_nts = set([('A', 'T'), ('T', 'A'), ('G', 'C'), ('C', 'G')])
opp_strand_dict = {'A':'T', 'G':'C', 'T':'A', 'C':'G'}

valid_nts = set(['A','T','C','G'])





def parse_parameters():
    """
    Parse the parameters into a dict, etc.
    """
    #    if len(sys.argv) == 1:
    #        print __doc__
    #        sys.exit(2)

    long_options_list = ['gf=',"gf_list=", 'vgf=', 'ssf=', 'out=', 'vbim=', 'N=', 'ssf_format=', 'gmdir=', 'indiv_list=',
                         'gf_format=', 'maf=', 'skip_coordination', 'skip_ambiguous', 'check_mafs', 'h', 'help',
                         'debug', 'keep_all', 'filter_sldsc=', "sigma2_trait=", "zscore", "FUNCT_separate",
                         "FUNCT_FILE="]

    p_dict = {'gf': None,"gf_list":None, 'vgf': None, 'ssf': None, 'out': None, 'vbim': None, 'N': None, 'ssf_format': 'STANDARD',
              'gmdir': None,
              'indiv_list': None, 'gf_format': 'PLINK', 'maf': 0.01, 'skip_coordination': False,
              'skip_ambiguous': False, 'debug': False, 'check_mafs': False, 'keep_all': False, 'filter_sldsc': 0,
              "sigma2_trait": 1, "zscore": False, "FUNCT_separate": False, "FUNCT_FILE": None}

    if len(sys.argv) == 1:
        print __doc__
    elif len(sys.argv) > 1:
        try:
            opts, args = getopt.getopt(sys.argv[1:], "h", long_options_list)

        except:
            print "Some problems with parameters.  Please read the usage documentation carefully."
            print "Use the -h option for usage information."
            #             traceback.print_exc()
            #             print __doc__
            sys.exit(2)

        for opt, arg in opts:
            if opt == "-h" or opt == "--h" or opt == '--help':
                print __doc__
                sys.exit(0)
            elif opt in ("--gf"):
                p_dict['gf'] = arg
            elif opt in ("--gf_list"):
                p_dict['gf_list'] = arg
            elif opt in ("--vgf"):
                p_dict['vgf'] = arg
            elif opt in ("--ssf"):
                p_dict['ssf'] = arg
            elif opt in ("--out"):
                p_dict['out'] = arg
            elif opt in ("--vbim"):
                p_dict['vbim'] = arg
            elif opt in ("--gmdir"):
                p_dict['gmdir'] = arg
            elif opt in ("--indiv_list"):
                p_dict['indiv_list'] = arg
            elif opt in ("--gf_format"):
                p_dict['gf_format'] = arg
            elif opt in ("--skip_coordination"):
                p_dict['skip_coordination'] = True
            elif opt in ("--skip_ambiguous"):
                p_dict['skip_ambiguous'] = True
            elif opt in ("--check_mafs"):
                p_dict['check_mafs'] = True
            elif opt in ("--keep_all"):
                p_dict['keep_all'] = True
            elif opt in ("--zscore"):
                p_dict['zscore'] = True
            elif opt in ("--FUNCT_separate"):
                p_dict['FUNCT_separate'] = True
            elif opt in ("--FUNCT_FILE"):
                p_dict['FUNCT_FILE'] = arg
            elif opt in ("--maf"):
                p_dict['maf'] = float(arg)
            elif opt in ("--filter_sldsc"):
                p_dict['filter_sldsc'] = float(arg)
            elif opt in ("--debug"):
                p_dict['debug'] = True
            elif opt in ("--ssf_format"):
                p_dict['ssf_format'] = arg
            elif opt in ("--N"):
                p_dict['N'] = int(arg)
            elif opt in ("--sigma2_trait"):
                p_dict['sigma2_trait'] = float(arg)
            else:
                print "Unkown option:", opt
                print "Use -h option for usage information."
                sys.exit(2)
    else:
        print __doc__
        sys.exit(0)
    return p_dict


def _get_chrom_dict_(loci, chromosomes):
    chr_dict = {}
    for chrom in chromosomes:
        chr_str = 'chrom_%d' % chrom
        chr_dict[chr_str] = {'sids': [], 'snp_indices': [], 'positions': [], 'nts': []}

    for i, l in enumerate(loci):
        chrom = l.chromosome
        pos = l.bp_position
        chr_str = 'chrom_%d' % chrom
        chr_dict[chr_str]['sids'].append(l.name)
        #         chr_dict[chr_str]['sids'].append('%d_%d'%(chrom,pos))
        chr_dict[chr_str]['snp_indices'].append(i)
        chr_dict[chr_str]['positions'].append(pos)
        chr_dict[chr_str]['nts'].append([l.allele1, l.allele2])

    print 'Bim file information loaded'
    return chr_dict



def _parse_plink_snps_(genotype_file, snp_indices):
    plinkf = plinkfile.PlinkFile(genotype_file)
    samples = plinkf.get_samples()
    num_individs = len(samples)
    num_snps = len(snp_indices)
    raw_snps = sp.empty((num_snps,num_individs),dtype='int8')
    #If these indices are not in order then we place them in the right place while parsing SNPs.
    snp_order = sp.argsort(snp_indices)
    ordered_snp_indices = list(snp_indices[snp_order])
    ordered_snp_indices.reverse()
    print 'Iterating over file to genotypes'
    snp_i = 0
    next_i = ordered_snp_indices.pop()
    line_i = 0
    max_i = ordered_snp_indices[0]
    while line_i <= max_i:
        if line_i < next_i:
            plinkf.next()
        elif line_i==next_i:
            line = plinkf.next()
            snp = sp.array(line, dtype='int8')
            bin_counts = line.allele_counts()
            if bin_counts[-1]>0:
                mode_v = sp.argmax(bin_counts[:2])
                snp[snp==3] = mode_v
            s_i = snp_order[snp_i]
            raw_snps[s_i]=snp
            if line_i < max_i:
                next_i = ordered_snp_indices.pop()
            snp_i+=1
        line_i +=1
    plinkf.close()
    assert snp_i==len(raw_snps), 'Failed to parse SNPs?'
    num_indivs = len(raw_snps[0])
    freqs = sp.sum(raw_snps,1, dtype='float32')/(2*float(num_indivs))
    return raw_snps, freqs


def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]


def parse_sum_stats_standard_ldscore(filename=None,
                                     bimfile_name=None,
                             hdf5_file_name=None,
                             n=None,outfile=None,filter=0,FUNCT_separate=True,FUNCT_FILE=False):
    """
    Summary statistics file must have the following columns with the following header (order is not important)
    CHR BP SNP A1 A2 BETA Z P
    """
    hdf5_file = h5py.File(hdf5_file_name, 'a')

    if bimfile_name is not None:
        bimf1 = re.sub(r"\[1:22\]", "[0-9]", bimfile_name)
        bimf2 = re.sub(r"\[1:22\]", "[0-2][0-9]", bimfile_name)

        print 'Parsing SNP list'
        valid_sids = set()
        bimfile_list=glob.glob(bimf1+".bim")+glob.glob(bimf2+".bim")
        bimfile_list.sort(key=natural_keys)
        for bimfile in bimfile_list:
            print 'Parsing bim file: %s' % bimfile
            with open(bimfile) as f:
                for line in f:
                    l = line.split()
                    valid_sids.add(l[1])
        print 'Total %d SNPs in all genotypes files.' % (len(valid_sids))

    chrom_dict = {}

    funct_dict = {}
    if FUNCT_separate:
        print 'Parsing the SNP-Heritability file: %s' % FUNCT_FILE
        with open(FUNCT_FILE) as f:
            print f.next()
            for line in f:
                l = (line.strip()).split()
                if float(l[1]) > 0:
                    if l[0] in valid_sids:
                        l = (line.strip()).split()
                        funct_dict[l[0]] = float(l[1])


    if FUNCT_separate:
        if len(valid_sids) == len(funct_dict.keys()):
            pass
        else:
            print("Get intersection between SNPs in genotypes file and functional enrichments file.")
            valid_sids = set(set(valid_sids) & set(funct_dict.keys()))
            print("%d found in both files"%len(valid_sids))

    snps_inf=[]
    ps_inf=[]
    beta_inf=[]
    beta_norm = []
    n_inf=0
    n_snps_tot=0
    n_snps_h2pos=0
    print 'Parsing the summary statistics file: %s' % filename
    with open(filename) as f:
        header = (f.next().strip()).split()
        try:
            idx_CHR= header.index("CHR")
            idx_BP = header.index("BP")
            idx_SNP = header.index("SNP")
            idx_BETA = header.index("BETA")
            idx_Z = header.index("Z")
            idx_P = header.index("P")
            idx_A1 = header.index("A1")
            idx_A2 = header.index("A2")

            for line in f:
                l = (line.strip()).split()
                try:
                    chrom = int(re.sub("[^0-9]","",l[idx_CHR]))
                    sid = l[idx_SNP]
                    if sid in valid_sids:
                        n_snps_tot+=1
                        if not chrom in chrom_dict.keys():
                            chrom_dict[chrom] = {'ps': [], 'log_odds': [],'betas': [], 'nts': [], 'sids': [], 'positions': [], 'ld_score':[]}
                        if FUNCT_separate:
                            h2snp_i = float(funct_dict[sid])
                        if h2snp_i > filter:
                            if sp.isinf(stats.norm.ppf(float(l[idx_P]) / 2.0) / sp.sqrt(n)):
                                n_inf += 1
                                snps_inf.append(sid)
                                ps_inf.append(float(l[idx_P]))
                                beta_inf.append(float(l[idx_BETA]))
                                chrom_dict[chrom]['sids'].append(sid)
                                chrom_dict[chrom]['positions'].append(l[idx_BP])
                                pval = float(l[idx_P])
                                chrom_dict[chrom]['ps'].append(pval)
                                nt = [l[idx_A1], l[idx_A2]]
                                chrom_dict[chrom]['nts'].append(nt)
                                raw_beta = float(l[idx_BETA])
                                chrom_dict[chrom]['log_odds'].append(raw_beta)
                                beta = sp.sign(float(l[idx_BETA]))*abs(float(l[idx_Z]))
                                beta_norm.append(beta / sp.sqrt(n))
                                chrom_dict[chrom]['betas'].append(beta / sp.sqrt(n))
                                chrom_dict[chrom]['ld_score'].append(float(funct_dict[sid]))
                            else:
                                n_snps_h2pos += 1
                                chrom_dict[chrom]['sids'].append(sid)
                                chrom_dict[chrom]['positions'].append(l[idx_BP])
                                pval = float(l[idx_P])
                                chrom_dict[chrom]['ps'].append(pval)
                                nt =  [l[idx_A1], l[idx_A2]]
                                chrom_dict[chrom]['nts'].append(nt)
                                raw_beta = float(l[idx_BETA])
                                chrom_dict[chrom]['log_odds'].append(raw_beta)
                                beta = sp.sign(raw_beta)*((-1)*stats.norm.ppf(pval / 2.0))
                                chrom_dict[chrom]['betas'].append(beta / sp.sqrt(n))
                                chrom_dict[chrom]['ld_score'].append(float(funct_dict[sid]))

                except ValueError:
                    print "Chromosome Not an integer!!! it is %s" % l[0][3:]

        except:
            print("Please check that the summary statistics contains a column for each of the follwing field:\n CHR (Chromosome)\n SNP \t SNP ID \n BP \t Physical position (base-pair) \n A1 \t Minor allele name (based on whole sample) \n A2 \t Major allele name \nP  \t Asymptotic p-value \n BETA \t Effect size \n Z \t Z-score \n")

    print '%d SNPs with p-value rounded to zero in summary statistics file' % n_inf
    if len(snps_inf)>0:
        print ' Writing SNPs with p-value rounded to zero in file: %s ' % outfile
        with open(outfile, 'w') as f:
            f.write("sid freqs pval beta \n")
            for sid,pval,beta,betnrm in it.izip(snps_inf,ps_inf,beta_inf,beta_norm):
                f.write("%s %0.4e %0.4e %0.4e\n"%(sid,pval,beta,betnrm))

    print 'SS file loaded, now sorting and storing in HDF5 file.'
    assert not 'sum_stats' in hdf5_file.keys(), 'Something is wrong with HDF5 file?'
    ssg = hdf5_file.create_group('sum_stats')
    num_snps = 0
    for chrom in chrom_dict.keys():
        print '%d SNPs on chromosome %d' % (len(chrom_dict[chrom]['sids']), chrom)
        sl = zip(chrom_dict[chrom]['sids'],chrom_dict[chrom]['positions'], chrom_dict[chrom]['nts'],
                 chrom_dict[chrom]['betas'], chrom_dict[chrom]['log_odds'],chrom_dict[chrom]['ld_score'],
                chrom_dict[chrom]['ps'])
        sl.sort()
        ps = []
        betas = []
        nts = []
        sids = []
        positions = []
        log_odds = []
        ld_score = []
        prev_pos = -1
        for sid, pos, nt, beta, lo, lsc,  p in sl:
            if pos == prev_pos:
                print 'duplicated position %d' % pos
                continue
            else:
                prev_pos = pos
            ps.append(p)
            betas.append(beta)
            nts.append(nt)
            sids.append(sid)
            positions.append(pos)
            log_odds.append(lo)
            ld_score.append(lsc)

        print 'Still %d SNPs on chromosome %d' % (len(ps), chrom)
        g = ssg.create_group('chrom_%d' % chrom)
        g.create_dataset('ps', data=sp.array(ps))
        g.create_dataset('betas', data=betas)
        g.create_dataset('log_odds', data=log_odds)
        g.create_dataset('ld_score', data=ld_score)
        num_snps += len(log_odds)
        g.create_dataset('nts', data=nts)
        g.create_dataset('sids', data=sids)
        g.create_dataset('positions', data=positions)
        hdf5_file.flush()
    hdf5_file.close()
    print '%d SNPs parsed from summary statistics file.' % num_snps




def coordinate_genot_ss(genotype_filename=None,
                        hdf5_file_name=None,
                        genetic_map_dir=None,
                        check_mafs=False,
                        min_maf=0.01,
                        skip_coordination=False, method=None, skip_ambiguous=False):
    """
    Assumes plink BED files.  Imputes missing genotypes.
    """
    hdf5_file = h5py.File(hdf5_file_name, 'a')
    ssf = hdf5_file['sum_stats']
    cord_data_g = hdf5_file.create_group('cord_data')

    bimf1 = re.sub(r"\[1:22\]", "[0-9]", genotype_filename)
    bimf2 = re.sub(r"\[1:22\]", "[0-2][0-9]", genotype_filename)
    bimfile_list = glob.glob(bimf1 + ".bim") + glob.glob(bimf2 + ".bim")
    bimfile_list.sort(key=natural_keys)

    print 'Coordinating summary statistics and genotypes'
    count_chr=0
    for bimfile in bimfile_list:
        count_chr+=1

        genotype_file = re.sub(r".bim", "", bimfile)

        if count_chr>1:
            hdf5_file = h5py.File(hdf5_file_name, 'a')
            cord_data_g = hdf5_file['cord_data']
            ssf = hdf5_file['sum_stats']

        #print("Processing file %s" % genotype_file)
        plinkf = plinkfile.PlinkFile(genotype_file)
        if count_chr==1:
            samples = plinkf.get_samples()
            num_individs = len(samples)
            #        num_individs = len(gf['chrom_1']['snps'][:, 0])
            #     Y = sp.array(gf['indivs']['phenotype'][...] == 'Case', dtype='int8')
            Y = [s.phenotype for s in samples]
            fids = [s.fid for s in samples]
            iids = [s.iid for s in samples]
            unique_phens = sp.unique(Y)
            if len(unique_phens) == 1:
                print 'Unable to find phenotype values.'
                has_phenotype = False
            elif len(unique_phens) == 2:
                cc_bins = sp.bincount(Y)
                assert len(cc_bins) == 2, 'Problems with loading phenotype'
                print 'Loaded %d controls and %d cases' % (cc_bins[0], cc_bins[1])
                has_phenotype = True
            else:
                print 'Found quantitative phenotype values'
                has_phenotype = True
            risk_scores = sp.zeros(num_individs)
            rb_risk_scores = sp.zeros(num_individs)
            if has_phenotype:
                hdf5_file.create_dataset('y', data=Y)

            hdf5_file.create_dataset('fids', data=fids)
            hdf5_file.create_dataset('iids', data=iids)


        num_common_snps = 0
        corr_list = []
        rb_corr_list = []



        # Figure out chromosomes and positions by looking at SNPs.
        loci = plinkf.get_loci()
        plinkf.close()
        gf_chromosomes = [l.chromosome for l in loci]

        chromosomes = sp.unique(gf_chromosomes)
        if len(chromosomes)>1:
            print( "Warning: There are more that one chromosomes in this plink file. This code assumes one file per chromosome.")
        #print( "Number of chromosomes in file %d" % len(chromosomes)) ### there should be only one
        chromosomes.sort()
        chrom=chromosomes
        #print chrom
        chr_str = 'chrom_%d' % chrom
        print 'Working on chromsome: %s' % chr_str

        chr_dict = _get_chrom_dict_(loci, chromosomes)

        tot_num_non_matching_nts = 0

        chrom_d = chr_dict[chr_str]
        try:
            ssg = ssf['chrom_%d' % chrom]
        except Exception, err_str:
            print err_str
            print 'Did not find chromsome in SS dataset.'
            print 'Continuing.'
            continue

        g_sids = chrom_d['sids']
        g_sid_set = set(g_sids)
        assert len(g_sid_set) == len(g_sids), 'Some duplicates?'
        ss_sids = ssg['sids'][...]
        ss_sid_set = set(ss_sids)
        assert len(ss_sid_set) == len(ss_sids), 'Some duplicates?'

        # Figure out filters:
        g_filter = sp.in1d(g_sids, ss_sids)
        ss_filter = sp.in1d(ss_sids, g_sids)

        # Order by SNP IDs
        g_order = sp.argsort(g_sids)
        ss_order = sp.argsort(ss_sids)

        g_indices = []
        for g_i in g_order:
            if g_filter[g_i]:
                g_indices.append(g_i)

        ss_indices = []
        for ss_i in ss_order:
            if ss_filter[ss_i]:
                ss_indices.append(ss_i)

        g_nts = chrom_d['nts']
        snp_indices = chrom_d['snp_indices']
        ss_nts = ssg['nts'][...]
        betas = ssg['betas'][...]
        log_odds = ssg['log_odds'][...]
        if method == 'STANDARD_FUNCT':
            ld_score = ssg['ld_score'][...]  ### S-LDSCORE
        #### Track allele flips indices ####
        ss_flips = sp.ones(len(ss_indices))
        assert not sp.any(sp.isnan(betas)), 'WTF?'
        # assert not sp.any(sp.isinf(betas)), 'WTF?'

        num_non_matching_nts = 0
        num_ambig_nts = 0
        ok_nts = []
        print 'Found %d SNPs present in both genotype and summary statistics datasets' % (len(g_indices))

        if 'freqs' in ssg.keys():
            ss_freqs = ssg['freqs'][...]
            ss_freqs_list = []

        ok_indices = {'g': [], 'ss': []}
        for g_i, ss_i in it.izip(g_indices, ss_indices):

            # Is the nucleotide ambiguous?
            # g_nt = [recode_dict[g_nts[g_i][0]],recode_dict[g_nts[g_i][1]]
            g_nt = [g_nts[g_i][0], g_nts[g_i][1]]

            if not skip_coordination:
                if not skip_ambiguous:
                    if tuple(g_nt) in ambig_nts:
                        num_ambig_nts += 1
                        tot_num_non_matching_nts += 1
                        continue

                if (not g_nt[0] in valid_nts) or (not g_nt[1] in valid_nts):
                    num_non_matching_nts += 1
                    tot_num_non_matching_nts += 1
                    continue

                ss_nt = ss_nts[ss_i]

                # Are the nucleotides the same?
                flip_nts = False
                os_g_nt = sp.array([opp_strand_dict[g_nt[0]], opp_strand_dict[g_nt[1]]])
                if not (sp.all(g_nt == ss_nt) or sp.all(os_g_nt == ss_nt)):
                    # Opposite strand nucleotides
                    flip_nts = (g_nt[1] == ss_nt[0] and g_nt[0] == ss_nt[1]) or (
                    os_g_nt[1] == ss_nt[0] and os_g_nt[0] == ss_nt[1])
                    if flip_nts:
                        betas[ss_i] = -betas[ss_i]
                        log_odds[ss_i] = -log_odds[ss_i]
                        ss_flips[ss_i] = -1
                        if 'freqs' in ssg.keys():
                            ss_freqs[ss_i] = 1 - ss_freqs[ss_i]
                    else:
                        #                     print "Nucleotides don't match after all?: g_sid=%s, ss_sid=%s, g_i=%d, ss_i=%d, g_nt=%s, ss_nt=%s" % \
                        #                         (g_sids[g_i], ss_sids[ss_i], g_i, ss_i, str(g_nt), str(ss_nt))
                        num_non_matching_nts += 1
                        tot_num_non_matching_nts += 1
                        continue

            # everything seems ok.
            ok_indices['g'].append(g_i)
            ok_indices['ss'].append(ss_i)
            ok_nts.append(g_nt)

        print '%d SNPs were excluded due to ambiguous nucleotides.' % num_ambig_nts
        print '%d SNPs were excluded due to non-matching nucleotides.' % num_non_matching_nts

        # Resorting by position
        positions = sp.array(chrom_d['positions'])[ok_indices['g']]
        order = sp.argsort(positions)
        ok_indices['g'] = list(sp.array(ok_indices['g'])[order])
        ok_indices['ss'] = list(sp.array(ok_indices['ss'])[order])
        positions = positions[order]

        # Parse SNPs
        snp_indices = sp.array(chrom_d['snp_indices'])
        snp_indices = snp_indices[ok_indices['g']]  # Pinpoint where the SNPs are in the file.
        raw_snps, freqs = _parse_plink_snps_(genotype_file, snp_indices)
        #print 'raw_snps.shape=', raw_snps.shape

        snp_stds = sp.sqrt(2 * freqs * (1 - freqs))  # sp.std(raw_snps, 1)
        snp_means = freqs * 2  # sp.mean(raw_snps, 1)

        betas = betas[ok_indices['ss']]
        log_odds = log_odds[ok_indices['ss']]
        if method == 'STANDARD_FUNCT':
            ld_score = ld_score[ok_indices['ss']]  #### S-LDSCORE
        ss_flips = ss_flips[ok_indices['ss']]  ### record flips
        ps = ssg['ps'][...][ok_indices['ss']]
        nts = sp.array(ok_nts)[order]
        sids = ssg['sids'][...][ok_indices['ss']]

        # Check SNP frequencies..
        if check_mafs and 'freqs' in ssg.keys():
            ss_freqs = ss_freqs[ok_indices['ss']]
            freq_discrepancy_snp = sp.absolute(ss_freqs - (1 - freqs)) > 0.15
            if sp.any(freq_discrepancy_snp):
                print 'Warning: %d SNPs appear to have high frequency discrepancy between summary statistics and validation sample' % sp.sum(
                    freq_discrepancy_snp)
                print freqs[freq_discrepancy_snp]
                print ss_freqs[freq_discrepancy_snp]

                # Filter freq_discrepancy_snps
                ok_freq_snps = sp.negative(freq_discrepancy_snp)
                raw_snps = raw_snps[ok_freq_snps]
                snp_stds = snp_stds[ok_freq_snps]
                snp_means = snp_means[ok_freq_snps]
                freqs = freqs[ok_freq_snps]
                ps = ps[ok_freq_snps]
                positions = positions[ok_freq_snps]
                nts = nts[ok_freq_snps]
                sids = sids[ok_freq_snps]
                betas = betas[ok_freq_snps]
                log_odds = log_odds[ok_freq_snps]
                if method == 'STANDARD_FUNCT':
                    ld_score = ld_score[ok_freq_snps]
                ss_flips = ss_flips[ok_freq_snps]

        # Filter minor allele frequency SNPs.
        maf_filter = (freqs > min_maf) * (freqs < (1 - min_maf))
        maf_filter_sum = sp.sum(maf_filter)
        n_snps = len(maf_filter)
        assert maf_filter_sum <= n_snps, "WTF?"
        if sp.sum(maf_filter) < n_snps:
            raw_snps = raw_snps[maf_filter]
            snp_stds = snp_stds[maf_filter]
            snp_means = snp_means[maf_filter]
            freqs = freqs[maf_filter]
            ps = ps[maf_filter]
            positions = positions[maf_filter]
            nts = nts[maf_filter]
            sids = sids[maf_filter]
            betas = betas[maf_filter]
            log_odds = log_odds[maf_filter]
            if method == 'STANDARD_FUNCT':
                ld_score = ld_score[maf_filter]

            print '%d SNPs with MAF < %0.3f were filtered' % (n_snps - maf_filter_sum, min_maf)

        print '%d SNPs were retained on chromosome %d.' % (maf_filter_sum, chrom)

        rb_prs = sp.dot(sp.transpose(raw_snps), log_odds)
        if has_phenotype:
            print 'Normalizing SNPs'
            snp_means.shape = (len(raw_snps), 1)
            snp_stds.shape = (len(raw_snps), 1)
            snps = (raw_snps - snp_means) / snp_stds
            assert snps.shape == raw_snps.shape, 'Aha!'
            snp_stds = snp_stds.flatten()
            snp_means = snp_means.flatten()
            prs = sp.dot(sp.transpose(snps), betas)
            corr = sp.corrcoef(Y, prs)[0, 1]
            corr_list.append(corr)
            print 'PRS correlation for chromosome %d was %0.4f' % (chrom, corr)
            rb_corr = sp.corrcoef(Y, rb_prs)[0, 1]
            rb_corr_list.append(rb_corr)
            print 'Raw effect sizes PRS correlation for chromosome %d was %0.4f' % (chrom, rb_corr)

        sid_set = set(sids)
        if genetic_map_dir is not None:
            genetic_map = []
            with gzip.open(genetic_map_dir + 'chr%d.interpolated_genetic_map.gz' % chrom) as f:
                for line in f:
                    l = line.split()
                    if l[0] in sid_set:
                        genetic_map.append(l[0])

        print 'Now storing coordinated data to HDF5 file from chrom_%d' % chrom[0]
        #print 'chrom_%d' % chrom[0]
        ofg = cord_data_g.create_group('chrom_%d' % chrom[0])
        ofg.create_dataset('raw_snps_ref', data=raw_snps, compression='lzf')
        ofg.create_dataset('snp_stds_ref', data=snp_stds)
        ofg.create_dataset('snp_means_ref', data=snp_means)
        ofg.create_dataset('freqs_ref', data=freqs)
        ofg.create_dataset('ps', data=ps)
        ofg.create_dataset('positions', data=positions)
        ofg.create_dataset('nts', data=nts)
        ofg.create_dataset('sids', data=sids)
        ofg.create_dataset('flips_ids', data=ss_flips)
        if genetic_map_dir is not None:
            ofg.create_dataset('genetic_map', data=genetic_map)
        # print 'Sum of squared effect sizes:', sp.sum(betas ** 2)
        #         print 'Sum of squared log odds:', sp.sum(log_odds ** 2)
        ofg.create_dataset('betas', data=betas)
        ofg.create_dataset('log_odds', data=log_odds)
        if method == 'STANDARD_FUNCT':
            ofg.create_dataset('ld_score', data=ld_score)
        ofg.create_dataset('log_odds_prs', data=rb_prs)
        if has_phenotype:
            risk_scores += prs
        rb_risk_scores += rb_prs
        num_common_snps += len(betas)
        #print "In chr %d there are %d betas and %d indices for flips" % (chrom, len(betas), len(ss_flips))
        hdf5_file.flush()
        hdf5_file.close()

        if has_phenotype:
            # Now calculate the prediction r^2
            corr = sp.corrcoef(Y, risk_scores)[0, 1]
            rb_corr = sp.corrcoef(Y, rb_risk_scores)[0, 1]
            print 'PRS R2 prediction accuracy for the whole genome was %0.4f (corr=%0.4f)' % (corr ** 2, corr)
            print 'Log-odds (effects) PRS R2 prediction accuracy for the whole genome was %0.4f (corr=%0.4f)' % (
            rb_corr ** 2, rb_corr)
    print 'There were %d SNPs in common' % num_common_snps
    print 'In all, %d SNPs were excluded due to nucleotide issues.' % tot_num_non_matching_nts
    print 'Done coordinating genotypes and summary statistics datasets.'


def main():
    p_dict = parse_parameters()
    print(p_dict)

    if p_dict['N'] is None:
        print 'Please specify an integer value for the sample size used to calculate the GWAS summary statistics.'

    print  'Preparing to parse summary statistics'
    if p_dict['vbim'] is not None:
        bimfile = p_dict['vbim']
    elif p_dict['vgf'] is not None:
        bimfile = p_dict['vgf']+'.bim'
    elif p_dict['gf'] is not None:
        bimfile = p_dict['gf']+'.bim'
    else:
        print 'Set of validation SNPs is missing!  Please specify either a validation PLINK genotype file, or a PLINK BIM file with the SNPs of interest.'
    if os.path.isfile(p_dict['out']):
        print 'Output file (%s) already exists!  Delete, rename it, or use a different output file.'%(p_dict['out'])
        raise Exception('Output file already exists!')

    # Open Hdf5 file
    #h5f = h5py.File(p_dict['out'],'a') ###
    parse_sum_stats_standard_ldscore(filename=p_dict['ssf'], bimfile_name=p_dict['gf'], hdf5_file_name=p_dict['out'], n=p_dict['N'],
                                    outfile=p_dict['out'] + "_snps_NaN.txt",FUNCT_FILE=p_dict["FUNCT_FILE"])


    coordinate_genot_ss(genotype_filename=p_dict['gf'], genetic_map_dir=p_dict['gmdir'], check_mafs=p_dict['check_mafs'],
                    hdf5_file_name=p_dict['out'], min_maf=p_dict['maf'], skip_coordination=p_dict['skip_coordination'],
                    method=p_dict['ssf_format'], skip_ambiguous=p_dict['skip_ambiguous'])

if __name__ == '__main__':
    main()

