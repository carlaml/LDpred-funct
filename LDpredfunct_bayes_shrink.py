#!/usr/bin/env python
import getopt
import sys
import traceback
import h5py
import scipy as sp
from scipy import linalg

import time

import itertools as it




def parse_parameters():
    """
    Parse the parameters into a dict, etc.
    """
#    if len(sys.argv) == 1:
#        print __doc__
#        sys.exit(2)

    long_options_list = ['coord=', 'ld_radius=', 'ld_prefix=', 'out=', 'N=', 'H2=',"alpha=",'ld_score=',"shrink=","M=","SVD=","method="]

    p_dict = {'coord':None, 'ld_radius':None, 'ld_prefix':None, 'out':None, 'N':None, 'H2':None,'alpha':None,'ld_score':None,'shrink':None,"M":None,'SVD':None,"method":None}

    if len(sys.argv) > 1:
        try:
            opts, args = getopt.getopt(sys.argv[1:], "h", long_options_list)

        except:
            print "Some problems with usage.  Please read the usage documentation carefully."
            traceback.print_exc()
            print __doc__
            sys.exit(2)

        for opt, arg in opts:
            if opt == "-h" or opt=="--h" or opt=='--help':
                print __doc__
                sys.exit(0)
            elif opt in ("--coord"): p_dict['coord'] = arg
            elif opt in ("--ld_radius"): p_dict['ld_radius'] = int(arg)
            elif opt in ("--ld_prefix"): p_dict['ld_prefix'] = arg
            elif opt in ("--out"): p_dict['out'] = arg
            elif opt in ("--N"): p_dict['N'] = int(arg)
            elif opt in ("--H2"): p_dict['H2'] = float(arg)
            elif opt in ("--SVD"):
                p_dict['SVD'] = bool(arg)
            elif opt == "--alpha":
                p_dict['alpha'] = float(arg)
            elif opt in ("--shrink"):
                p_dict['shrink'] = float(arg)
            elif opt == "--ld_score":
                p_dict['ld_score'] = arg
            elif opt in ("--M"):
                p_dict['M'] = int(arg)
            elif opt in ("--method"):
                p_dict['method'] = arg
            else:
                print "Unkown option:", opt
                print "Use -h option for usage information."
                sys.exit(2)
    else:
        print __doc__
        sys.exit(0)
    return p_dict






def ldpred_inf_M_solve(beta_hats, h2=0.1, n=1000, inf_shrink_matrices=None,
               reference_ld_mats=None, genotypes=None, ld_window_size=100, verbose=False,m=None,SVD=False):
    """
    ## Since the input is only per chromosome, need to input externally the total number of SNPS
    Apply the infinitesimal shrink w LD (which requires LD information).

    If reference_ld_mats are supplied, it uses those, otherwise it uses the LD in the genotype data.

    If genotypes are supplied, then it assumes that beta_hats and the genotypes are synchronized.

    """
    if verbose:
        print 'Doing LD correction'
    t0 = time.time()
    num_betas = len(beta_hats)
    updated_betas = sp.empty(num_betas)

    for i, wi in enumerate(range(0, num_betas, ld_window_size)):
        start_i = wi
        stop_i = min(num_betas, wi + ld_window_size)
        curr_window_size = stop_i - start_i
        X = genotypes[start_i: stop_i]
        num_indivs = X.shape[1]
        D = sp.dot(X, X.T) / num_indivs

        A = ((m / h2) * sp.eye(curr_window_size) + (n / (1)) * D)
        updated_betas[start_i: stop_i] =linalg.solve(A,n*beta_hats[start_i: stop_i])

        if verbose:
            sys.stdout.write('\b\b\b\b\b\b\b%0.2f%%' % (100.0 * (min(1, float(wi + 1) / num_betas))))
            sys.stdout.flush()

    t1 = time.time()
    t = (t1 - t0)
    if verbose:
        print '\nIt took %d minutes and %0.2f seconds to perform the Infinitesimal LD shrink' % (t / 60, t % 60)
    return updated_betas



def ldpredfunct_solve_ldblocks(beta_hats, snp_h2, h2=0.1,  n=1000, genotypes=None, ld_window_size=100, verbose=False, Cvalue= None):
    """
    Apply the infinitesimal shrink w LD (which requires LD information).

    If reference_ld_mats are supplied, it uses those, otherwise it uses the LD in the genotype data.

    If genotypes are supplied, then it assumes that beta_hats and the genotypes are synchronized.

    """
    #print("Shrink constant c: "+str(Cvalue))
    if verbose:
        print 'Doing LD correction'
    t0 = time.time()
    num_betas = len(beta_hats)
    updated_betas = sp.empty(num_betas)

    for i, wi in enumerate(range(0, num_betas, ld_window_size)):
        start_i = wi
        stop_i = min(num_betas, wi + ld_window_size)
        curr_window_size = stop_i - start_i
        X = genotypes[start_i: stop_i]
        num_indivs = X.shape[1]
        D = sp.dot(X, X.T) / num_indivs


        Winv = (1 / Cvalue) * (snp_h2[start_i: stop_i]** (- 1) * sp.eye(curr_window_size))
        A = (Winv + (n / (1)) * D)
        updated_betas[start_i: stop_i] =linalg.solve(A,n*beta_hats[start_i: stop_i])
        if verbose:
            sys.stdout.write('\b\b\b\b\b\b\b%0.2f%%' % (100.0 * (min(1, float(wi + 1) / num_betas))))
            sys.stdout.flush()
    t1 = time.time()
    t = (t1 - t0)
    if verbose:
        print '\nIt took %d minutes and %0.2f seconds to perform the Infinitesimal LD shrink' % (t / 60, t % 60)
    return updated_betas


def get_LDpred_single_ld_table(wi,snps, ld_window_size=200):
    """
    Calculates LD tables, and the LD score in one go...
    """
    m, n = snps.shape
    start_i = wi
    stop_i = min(m, wi + ld_window_size)
    X = snps[start_i: stop_i]
    D = sp.dot(X, X.T) / n
    return D


def ldpred_funct_genomewide(data_file=None, ld_radius = None, out_file_prefix=None,
                          n=None, h2=None, Cval=None,method="ldscore_solve",m=None,SVD=False):
    """
    Calculate LDpred for a genome
    """

    df = h5py.File(data_file,'r')
    has_phenotypes=False
    if 'y' in df.keys():
        'Validation phenotypes found.'
        y = df['y'][...]  # Phenotype
        num_individs = len(y)
        risk_scores_pval_derived = sp.zeros(num_individs)
        has_phenotypes=True


    results_dict = {}
    num_snps = 0
    sum_beta2s = 0
    cord_data_g = df['cord_data']

    tot_snp_sum = 0 #### <<-----
    totsum_freqs = 0  #### <<-----

    for chrom_str in cord_data_g.keys():
        g = cord_data_g[chrom_str]
        #print g.keys()
        betas = g['betas'][...]
        n_snps = len(betas)
        num_snps += n_snps
        sum_beta2s += sp.sum(betas ** 2)
        snp_stds = g['snp_stds_ref'][...]  #### <<-----
        snp_stds = snp_stds.flatten() #### <<-----
        ok_snps_filter = snp_stds > 0 #### <<-----
        # snp_freqs = g['freqs_ref'][...]  #### <<-----
        # snp_freqs = snp_freqs[ok_snps_filter]  ### <<----
        #print ( "In chromosome %s, we have %d SNPs"%(chrom_str,n_snps))
        if method == "ldscore_solve" or method == "ldscore+ldpredinf_solve":
            snp_h2=g['ld_score'][...]
            snp_h2 = snp_h2[ok_snps_filter]  ### <<----
            tot_snp_sum += sp.sum(snp_h2)

    if m==None:
        m=num_snps

    print ( "Total number of SNPs: %d"%m)

    if ld_radius==None:
        ld_window=round((0.15/100)*m)- round((0.15/100)*m)%100 ### set window size to a close number
    else:
        ld_window=2 * ld_radius

    print ("LD window size: %d" % ld_window)

    if method == "ldscore_solve" or method == "ldscore+ldpredinf_solve":
        if Cval==None:
            try:
                Cval = h2 / tot_snp_sum
            except ValueError:
                print ("Error computing normalizing constant: h2g is %0.4e and sum of positive enrichments across all SNPs is %0.4e" %(h2, tot_snp_sum) )
        print("Normalizing constant for an h2g of %0.4e is set to be %0.12e"%(h2,Cval))

    chi_square_lambda = sp.mean(n * sum_beta2s / float(num_snps))
    assert chi_square_lambda>1, 'Something is wrong with the GWAS summary statistics, parsing of them, or the given GWAS sample size (N). Lambda (the mean Chi-square statistic) is too small.  '


    if out_file_prefix:
        #Preparing output files
        raw_effect_sizes = []
        ldpred_effect_sizes = []
        ldpred_effect_sizes_ldsc = []
        sids = []
        chromosomes = []
        positions = []
        nts = []

    for chrom_str in cord_data_g.keys():
        g = cord_data_g[chrom_str]
        # if has_phenotypes:
        #     if 'raw_snps_val' in g.keys():
        #         raw_snps = g['raw_snps_val'][...]
        #     else:
        #         raw_snps = g['raw_snps_ref'][...]

        raw_snps = g['raw_snps_ref'][...]
        snp_stds = g['snp_stds_ref'][...]
        snp_means = g['snp_means_ref'][...]

        # Filter monomorphic SNPs
        snp_stds = snp_stds.flatten()
        ok_snps_filter = snp_stds > 0
        snp_stds = snp_stds[ok_snps_filter]
        pval_derived_betas = g['betas'][...]
        pval_derived_betas = pval_derived_betas[ok_snps_filter]

        #snp_freqs = g['freqs_ref'][...]  #### <<-----
        #snp_freqs = snp_freqs[ok_snps_filter]  ### <<----

        n_snps = len(raw_snps)
        snp_means.shape = (n_snps, 1)
        snp_stds.shape = (n_snps, 1)

        snps = sp.array((raw_snps - snp_means) / snp_stds, dtype='float32')

        if method=="ldscore_solve" or method == "ldscore+ldpredinf_solve":
            snp_h2=g['ld_score'][...]
            snp_h2 = snp_h2[ok_snps_filter]  ### <<----


        if out_file_prefix:
            chromosomes.extend([chrom_str]*len(pval_derived_betas))
            positions.extend(g['positions'][...])
            sids.extend(g['sids'][...])
            raw_effect_sizes.extend(g['log_odds'][...])
            nts.extend(g['nts'][...])

        n_snps = len(pval_derived_betas)

        if h2 is not None:
            h2_chrom = h2 * (n_snps / float(num_snps))
            #print h2_chrom ### Note that here h2_chrom is h2


        print 'Calculating posterior means for chromosome %s'%((chrom_str.split('_'))[1])

        if method == "ldpredinf_solve":
            updated_betas = ldpred_inf_M_solve(pval_derived_betas, genotypes=snps,
                                         h2=h2, n=n, ld_window_size=ld_window, verbose=True, m=m, SVD=SVD)
        elif method == "ldscore_solve":
            updated_betas = ldpredfunct_solve_ldblocks(pval_derived_betas, snp_h2, genotypes=snps,
                                               h2=h2_chrom, n=n, ld_window_size=ld_window, verbose=False,
                                               Cvalue=Cval)

        updated_betas = updated_betas / (snp_stds.flatten())
        ldpred_effect_sizes.extend(updated_betas)


    if method == "ldscore_solve":
        weights_out_file = '%s_LDpred-inf-ldscore.txt' % (out_file_prefix)
        print( " Writing file %s" % weights_out_file)
        with open(weights_out_file, 'w') as f:
            f.write('chrom    pos    sid    nt1    nt2    raw_beta    ldpred_inf_beta\n')
            for chrom, pos, sid, nt, raw_beta, ldpred_beta in it.izip(chromosomes, positions, sids, nts,
                                                                      raw_effect_sizes, ldpred_effect_sizes):
                nt1, nt2 = nt[0], nt[1]
                f.write('%s    %d    %s    %s    %s    %0.4e    %0.4e\n' % (
                chrom, pos, sid, nt1, nt2, raw_beta, ldpred_beta))
    elif method == "ldpredinf_solve":
        weights_out_file = '%s_LDpred-inf.txt' % (out_file_prefix)
        print( " Writing file %s" % weights_out_file)

        with open(weights_out_file, 'w') as f:
            f.write('chrom    pos    sid    nt1    nt2    raw_beta    ldpred_inf_beta\n')
            for chrom, pos, sid, nt, raw_beta, ldpred_beta in it.izip(chromosomes, positions, sids, nts,
                                                                      raw_effect_sizes, ldpred_effect_sizes):
                nt1, nt2 = nt[0], nt[1]
                f.write('%s    %d    %s    %s    %s    %0.4e    %0.4e\n' % (
                chrom, pos, sid, nt1, nt2, raw_beta, ldpred_beta))





def main():
    p_dict = parse_parameters()
    print p_dict


    print """
Note: For maximal accuracy all SNPs with LDpred weights should be included in the validation data set.
If they are a subset of the validation data set, then we suggest recalculate LDpred for the overlapping SNPs.
"""

    ldpred_funct_genomewide(data_file=p_dict['coord'], out_file_prefix=p_dict['out'], ld_radius=p_dict['ld_radius'],
                          n=p_dict['N'],h2=p_dict['H2'],Cval=p_dict['shrink'],m=p_dict['M'],method=p_dict["method"])


if __name__ == '__main__':
    main()


