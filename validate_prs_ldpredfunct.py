#!/usr/bin/env python

import getopt
import sys
import os
import traceback
import scipy as sp
from scipy import linalg
from plinkio import plinkfile
import itertools as it
import time
import h5py
import math
import re
import glob

ok_nts = ['A','T','C','G']
opp_strand_dict = {'A':'T', 'G':'C', 'T':'A', 'C':'G'}
non_zero_chromosomes =set()





def parse_parameters():
    """
    Parse the parameters into a dict, etc.
    """
    long_options_list = ['vgf=', 'rf=', 'res_format=', 'out=', 'indiv_filter=', 'split_by_chrom', 'pf=', 'pf_format=', 'cov_file=',
                         'pcs_file=', 'PS=', 'TS=',"method=", 'adjust_for_sex','adjust_for_covariates', 'adjust_for_pcs','flips', 'h','help','Q_tresh=',"beta_thresh=","beta_thresh_file="]

    p_dict = {'vgf':None, 'rf':None, 'out':None, 'res_format':'LDPRED', 'indiv_filter':None, 'split_by_chrom':False,
              'pf':None, 'pf_format':'STANDARD', 'cov_file':None, "method":None,'pcs_file':None, 'PS':[1,0.3,0.1,0.03,0.01,0.003,0.001],
              'TS':[1,0.3,0.1,0.03,0.01,0.003,0.001,3*1E-4,1E-4,3*1E-5,1E-5,1E-6,1E-7,1E-8],
              'adjust_for_sex':False, 'adjust_for_covariates':False, 'adjust_for_pcs':False,'Q_tresh':[1.7,2.5,4.3],'flips':False,"beta_thresh":False,"beta_thresh_file":None}

    if len(sys.argv) > 1:
        try:
            opts, args = getopt.getopt(sys.argv[1:], "h", long_options_list)

        except:
            print "Some problems with parameters.  Please read the usage documentation carefully."
            print "Use the -h option for usage information."
#             traceback.print_exc()
#             print __doc__
            sys.exit(2)

        for opt, arg in opts:
            if opt == "-h" or opt=="--h" or opt=='--help':
                print __doc__
                sys.exit(0)
            elif opt in ("--vgf"): p_dict['vgf'] = arg
            elif opt in ("--rf"): p_dict['rf'] = arg
            elif opt in ("--res_format"): p_dict['res_format'] = arg
            elif opt in ("--indiv_filter"): p_dict['indiv_filter'] = arg
            elif opt in ("--out"): p_dict['out'] = arg
            elif opt in ("--method"): p_dict['method'] = arg

            elif opt in ("--split_by_chrom"): p_dict['split_by_chrom'] = True
            elif opt in ("--PS"): p_dict['PS'] = map(float,arg.split(','))
            elif opt in ("--beta_thresh"): p_dict['beta_thresh'] = map(float,arg.split(','))
            elif opt in ("--beta_thresh_file"): p_dict['beta_thresh_file'] =arg

            elif opt in ("--TS"): p_dict['TS'] = map(float,arg.split(','))
            elif opt in ("--Q_tresh"):
                p_dict['Q_tresh'] = map(float, arg.split(','))
            elif opt in ("--pf"): p_dict['pf'] = arg
            elif opt in ("--pf_format"): p_dict['pf_format'] = arg
            elif opt in ("--cov_file"): p_dict['cov_file'] = arg
            elif opt in ("--pcs_file"): p_dict['pcs_file'] = arg
            elif opt in ("--flips"): p_dict['flips'] = True
            elif opt in ("--adjust_for_sex"): p_dict['adjust_for_sex'] = True
            elif opt in ("--adjust_for_covariates"): p_dict['adjust_for_covariates'] = True
            elif opt in ("--adjust_for_pcs"): p_dict['adjust_for_pcs'] = True
            else:
                print "Unkown option:", opt
                print "Use -h option for usage information."
                sys.exit(2)
    else:
        print __doc__
        sys.exit(0)
    return p_dict



def parse_ldpred_res_bins(file_name,K_bins=1):
    rs_id_map = {}
    """
    chrom    pos    sid    nt1    nt2    raw_beta    ldpred_inf_beta    ldpred_beta
    1, 798959, rs11240777, C, T, -1.1901e-02, 3.2443e-03, 2.6821e-04
    """
    upd_pval_beta_sq_list=[]
    sum_pval_beta_sq=0
    with open(file_name,'r') as f:
        f.next()
        for line in f:
            l = line.split()
            chrom_str = l[0]
            chrom = int(chrom_str[6:])
            pos = int(l[1])
            rs_id = l[2].strip()
            nt1 = l[3].strip()
            nt2 = l[4].strip()
            nts = [nt1,nt2]
            raw_beta = float(l[5])
            upd_pval_beta = float(l[6])
            upd_pval_beta_sq_list.append(float(l[6])**2)
            sum_pval_beta_sq+=float(l[6])**2
            rs_id_map[rs_id] = {'chrom':chrom, 'pos':pos, 'nts':nts, 'raw_beta':raw_beta,
                                'upd_pval_beta':upd_pval_beta}
    S_k=sum_pval_beta_sq/K_bins
    if K_bins>1:
        upd_pval_beta_sq_list.sort()
        upd_pval_beta_sq_list=sp.array(upd_pval_beta_sq_list)
        upd_pval_beta_sq_list.shape=(len(upd_pval_beta_sq_list),)
        cumsum_upd_pval_beta=sp.cumsum(upd_pval_beta_sq_list)


        S_k_list = [0]
        for k in range(1,K_bins):
            bk=S_k * k
            eff_bin =upd_pval_beta_sq_list[sp.where(cumsum_upd_pval_beta >= bk)]
            if min(eff_bin) in S_k_list: pass
            else:
                S_k_list.append(min(eff_bin))
        S_k_list.append(upd_pval_beta_sq_list[-1])
        S_k_list=sp.array(S_k_list)
        S_k_list.shape=(len(S_k_list),)
        rs_id_map["bins_extremes"]=S_k_list


    return rs_id_map


def get_prs_bins(genotype_file, rs_id_map, K_bins=1,phen_map=None,lasso=False,sets=False,verbose=False):
    plinkf = plinkfile.PlinkFile(genotype_file)
    samples = plinkf.get_samples()

    #1. Figure out indiv filter and get true phenotypes
    indiv_filter=sp.zeros(len(samples),dtype='bool8')
    true_phens = []
    iids = []
    if phen_map is not None:
        pcs = []
        sex = []
        covariates = []
        phen_iids = set(phen_map.keys())
        for samp_i, sample in enumerate(samples):
            if sample.iid in phen_iids:
                indiv_filter[samp_i] = True
                true_phens.append(phen_map[sample.iid]['phen'])
                iids.append(sample.iid)
                if 'pcs' in phen_map[sample.iid].keys():
                    pcs.append(phen_map[sample.iid]['pcs'])
                if 'sex' in phen_map[sample.iid].keys():
                    sex.append(phen_map[sample.iid]['sex'])
                if 'covariates' in phen_map[sample.iid].keys():
                    #Temp hack...
#                     if phen_map[sample.iid]['sex']==1:
#                         covariates.append([phen_map[sample.iid]['covariates'][0],0])
#                     else:
#                         covariates.append([0,phen_map[sample.iid]['covariates'][0]])
                    covariates.append(phen_map[sample.iid]['covariates'])
        if len(pcs)>0:
            assert len(pcs)==len(true_phens), 'PC information missing for some individuals with phenotypes'
        if len(sex)>0:
            assert len(sex)==len(true_phens), 'Sex information missing for some individuals with phenotypes'
        if len(covariates)>0:
            assert len(covariates)==len(true_phens), 'Covariates missing for some individuals with phenotypes'
    else:
        for samp_i, sample in enumerate(samples):
            if sample.affection!=2:
                indiv_filter[samp_i] = True
                true_phens.append(sample.affection)
                iids.append(sample.iid)

    num_individs = sp.sum(indiv_filter)
    assert num_individs>0, 'Issues in parsing the phenotypes and/or PCs?'

    assert not sp.any(sp.isnan(true_phens)),'Phenotypes appear to have some NaNs, or parsing failed.'

    print '%d individuals have phenotype and genotype information.'%num_individs

    num_non_matching_nts = 0
    num_flipped_nts = 0

    raw_effects_prs = sp.zeros(num_individs)
    pval_derived_effects_prs = sp.zeros(num_individs)
    pval_derived_effects_prs_lasso = sp.zeros(num_individs)

    bins_prs_dict={}
    if K_bins>1:
        bk=1
        while bk <= K_bins:
            bins_prs_dict["prs_bin_%d" % bk] = sp.zeros(num_individs)
            bk += 1
        

    #Sets
    pval_derived_effects_prs_high = sp.zeros(num_individs)
    pval_derived_effects_prs_lasso_high = sp.zeros(num_individs)
    pval_derived_effects_prs_low = sp.zeros(num_individs)
    pval_derived_effects_prs_lasso_low = sp.zeros(num_individs)
    #If these indices are not in order then we place them in the right place while parsing SNPs.
    print 'Iterating over BED file to calculate risk scores.'
    locus_list = plinkf.get_loci()
    snp_i = 0
    if K_bins > 1:
        bins_bounds=rs_id_map["bins_extremes"]
    #print bins_bounds
    for locus, row in it.izip( locus_list, plinkf):
        upd_pval_beta = 0
        try:
            #Check rs-ID
#             sid = '%d_%d'%(locus.chromosome,locus.bp_position)
            sid = locus.name
            rs_info = rs_id_map[sid]
        except Exception: #Move on if rsID not found.
            continue

        if rs_info['upd_pval_beta']==0:
            continue

        #Check whether the nucleotides are OK, and potentially flip it.
        ss_nt = rs_info['nts']
        g_nt =  [locus.allele1,locus.allele2]
        flip_nts = False
        os_g_nt = sp.array([opp_strand_dict[g_nt[0]], opp_strand_dict[g_nt[1]]])
        if not (sp.all(g_nt == ss_nt) or sp.all(os_g_nt == ss_nt)):
            # Opposite strand nucleotides
            flip_nts = (g_nt[1] == ss_nt[0] and g_nt[0] == ss_nt[1]) or (os_g_nt[1] == ss_nt[0] and os_g_nt[0] == ss_nt[1])
            if flip_nts:
                raw_beta = -rs_info['raw_beta']
                upd_pval_beta = -rs_info['upd_pval_beta']
                num_flipped_nts+=1
                if lasso:
                    upd_pval_beta_lasso = -rs_info['upd_pval_beta_lasso']
                    if sets:
                        upd_pval_beta_high = -rs_info['upd_pval_beta_high']
                        upd_pval_beta_lasso_high = -rs_info['upd_pval_beta_lasso_high']
                        upd_pval_beta_low = -rs_info['upd_pval_beta_low']
                        upd_pval_beta_lasso_low = -rs_info['upd_pval_beta_lasso_low']
            else:
                #print "Nucleotides don't match after all?: sid=%s, g_nt=%s, ss_nt=%s" % (locus.name, str(g_nt), str(ss_nt))
                num_non_matching_nts += 1
                continue
        else:
            raw_beta = rs_info['raw_beta']
            upd_pval_beta = rs_info['upd_pval_beta']
            if lasso:
                upd_pval_beta_lasso = rs_info['upd_pval_beta_lasso']
                if sets:
                    upd_pval_beta_high = rs_info['upd_pval_beta_high']
                    upd_pval_beta_lasso_high = rs_info['upd_pval_beta_lasso_high']
                    upd_pval_beta_low = rs_info['upd_pval_beta_low']
                    upd_pval_beta_lasso_low = rs_info['upd_pval_beta_lasso_low']

        #Parse SNP, and fill in the blanks if necessary.
        snp = sp.array(row, dtype='int8')[indiv_filter]
        bin_counts = row.allele_counts()
        if bin_counts[-1]>0:
            mode_v = sp.argmax(bin_counts[:2])
            snp[snp==3] = mode_v

        #Normalize SNP
#         n_snp = (snp - sp.mean(snp))/sp.std(snp)
#         print(upd_pval_beta**2)
#         print sp.where(bins_bounds>=upd_pval_beta**2)
#         print sp.where(bins_bounds>=upd_pval_beta**2)[0][0]

        #Update scores and move on.
        raw_effects_prs += snp*raw_beta
        assert not sp.any(sp.isnan(raw_effects_prs)),'Raw effects PRS is corrupted'
        snpi_b = snp * upd_pval_beta
        pval_derived_effects_prs += snpi_b
        if K_bins > 1:
            bin_number = sp.where(bins_bounds >= upd_pval_beta ** 2)[0][0]
            bins_prs_dict["prs_bin_%d" % bin_number] += snpi_b
        assert not sp.any(sp.isnan(pval_derived_effects_prs)),'Weighted effects PRS is corrupted'

        if verbose:

            if snp_i>0 and snp_i%500000==0:
                print("PRS using %d SNPS"%snp_i)
                #print 'Number of non-matching NTs: %d'%num_non_matching_nts
                raw_eff_r2 = (sp.corrcoef(raw_effects_prs, true_phens)[0,1])**2
                pval_eff_r2  = (sp.corrcoef(pval_derived_effects_prs, true_phens)[0,1])**2
                print 'Raw effects PRS r2: %0.4f'%raw_eff_r2
                print 'Weigted effects PRS r2: %0.4f'%pval_eff_r2
                if lasso:
                    pval_eff_r2_lasso = (sp.corrcoef(pval_derived_effects_prs_lasso, true_phens)[0, 1]) ** 2
                    print 'Weigted effects PRS Lasso r2: %0.4f' % pval_eff_r2_lasso

                    if sets:
                        pval_eff_r2_high = (sp.corrcoef(pval_derived_effects_prs_high, true_phens)[0, 1]) ** 2
                        print 'Weigted effects HIGH PRS r2: %0.4f' % pval_eff_r2_high
                        pval_eff_r2_lasso_high = (sp.corrcoef(pval_derived_effects_prs_lasso_high, true_phens)[0, 1]) ** 2
                        print 'Weigted effects HIGH PRS Lasso r2: %0.4f' % pval_eff_r2_lasso_high

                        pval_eff_r2_low = (sp.corrcoef(pval_derived_effects_prs_low, true_phens)[0, 1]) ** 2
                        print 'Weigted effects LOW PRS r2: %0.4f' % pval_eff_r2_low
                        pval_eff_r2_lasso_low = (sp.corrcoef(pval_derived_effects_prs_lasso_low, true_phens)[0, 1]) ** 2
                        print 'Weigted effects LOW PRS Lasso r2: %0.4f' % pval_eff_r2_lasso_low

        snp_i +=1

    plinkf.close()


    print "DONE!"
    print 'Number of non-matching NTs: %d'%num_non_matching_nts
    print 'Number of flipped NTs: %d'%num_flipped_nts
    raw_eff_corr = sp.corrcoef(raw_effects_prs, true_phens)[0,1]
    raw_eff_r2 = raw_eff_corr**2
    pval_eff_corr = sp.corrcoef(pval_derived_effects_prs, true_phens)[0,1]
    pval_eff_r2  = pval_eff_corr**2

    print 'Raw effects PRS correlation: %0.4f'%raw_eff_corr
    print 'Raw effects PRS r2: %0.4f'%raw_eff_r2
    print 'Weigted effects PRS correlation: %0.4f'%pval_eff_corr
    print 'Weigted effects PRS r2: %0.4f'%pval_eff_r2

    if lasso:
        pval_eff_corr_lasso = sp.corrcoef(pval_derived_effects_prs_lasso, true_phens)[0, 1]
        pval_eff_r2_lasso = pval_eff_corr_lasso ** 2
        print 'Weigted effects LASSO PRS correlation: %0.4f' % pval_eff_corr_lasso
        print 'Weigted effects LASSO PRS r2: %0.4f' % pval_eff_r2_lasso
        if sets:
            pval_eff_corr_high = sp.corrcoef(pval_derived_effects_prs_high, true_phens)[0, 1]
            pval_eff_r2_high = pval_eff_corr_high ** 2
            print 'Weigted effects HIGH PRS correlation: %0.4f' % pval_eff_corr_high
            print 'Weigted effects HIGH PRS r2: %0.4f' % pval_eff_r2_high
            pval_eff_corr_lasso_high = sp.corrcoef(pval_derived_effects_prs_lasso_high, true_phens)[0, 1]
            pval_eff_r2_lasso_high = pval_eff_corr_lasso_high ** 2
            print 'Weigted effects HIGH LASSO PRS correlation: %0.4f' % pval_eff_corr_lasso_high
            print 'Weigted effects HIGH LASSO PRS r2: %0.4f' % pval_eff_r2_lasso_high

            pval_eff_corr_low = sp.corrcoef(pval_derived_effects_prs_low, true_phens)[0, 1]
            pval_eff_r2_low = pval_eff_corr_low ** 2
            print 'Weigted effects LOW PRS correlation: %0.4f' % pval_eff_corr_low
            print 'Weigted effects LOW PRS r2: %0.4f' % pval_eff_r2_low
            pval_eff_corr_lasso_low = sp.corrcoef(pval_derived_effects_prs_lasso_low, true_phens)[0, 1]
            pval_eff_r2_lasso_low = pval_eff_corr_lasso_low ** 2
            print 'Weigted effects LOW LASSO PRS correlation: %0.4f' % pval_eff_corr_lasso_low
            print 'Weigted effects LOW LASSO PRS r2: %0.4f' % pval_eff_r2_lasso_low


    ret_dict = {'raw_effects_prs':raw_effects_prs.copy(), 'pval_derived_effects_prs':pval_derived_effects_prs.copy(),
                'true_phens':true_phens[:], 'iids':iids}

    if K_bins>1:
        bk=1
        while bk <= K_bins:
            ret_dict["pval_derived_effects_prs_bin_%d" % bk] = bins_prs_dict["prs_bin_%d" % bk].copy()
            bk += 1


    if len(pcs)>0:
        ret_dict['pcs'] = pcs
    if len(sex)>0:
        ret_dict['sex'] = sex
    if len(covariates)>0:
        ret_dict['covariates'] = covariates

    return ret_dict




def parse_phen_file(pf, pf_format):
    print pf
    phen_map ={}
    if pf!=None:
        if pf_format=='FAM':
            """
            Individual's family ID ('FID')
            Individual's within-family ID ('IID'; cannot be '0')
            Within-family ID of father ('0' if father isn't in dataset)
            Within-family ID of mother ('0' if mother isn't in dataset)
            Sex code ('1' = male, '2' = female, '0' = unknown)
            Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)
            """
            print 'Parsing phenotypes'
            with open(pf,'r') as f:
                for line in f:
                    l = line.split()
                    iid = l[1]
                    #iid = iid[:-4]
                    sex = int(l[4])
                    phen = float(l[5])
                    if sex!=0 and phen !=-9:
                        phen_map[iid] = {'phen':phen, 'sex':sex}

            iids = set(phen_map.keys())

        if pf_format=='STANDARD':
            """
            IID   PHE
            """
            print 'Parsing phenotypes'
            with open(pf,'r') as f:
                print(f)
                for line in f:
                    l = line.split()
                    iid = l[0]
                    phen = float(l[1])
                    phen_map[iid] = {'phen':phen}
            iids = set(phen_map.keys())


        elif pf_format=='S2':
            """
            IID Age Sex Height_Inches
            """
            with open(pf,'r') as f:
                print f.next()
                for line in f:
                    l = line.split()
                    iid = l[0]
                    age = float(l[1])
                    if l[2]=='Male':
                        sex = 1
                    elif l[2]=='Female':
                        sex = 2
                    else:
                        raise Exception('Sex missing')
                    phen = float(l[3])
                    phen_map[iid] = {'phen':phen, 'age':age, 'sex':sex}

    return phen_map

def parse_ldpred_res_beta2_up(file_name,betathresh=None):
    rs_id_map = {}
    """
    chrom    pos    sid    nt1    nt2    raw_beta    ldpred_inf_beta    ldpred_beta
    1, 798959, rs11240777, C, T, -1.1901e-02, 3.2443e-03, 2.6821e-04
    """
    nsnps=0
    with open(file_name,'r') as f:
        f.next()
        for line in f:
            l = line.split()
            chrom_str = l[0]
            chrom = int(chrom_str[6:])
            pos = int(l[1])
            rs_id = l[2].strip()
            nt1 = l[3].strip()
            nt2 = l[4].strip()
            nts = [nt1,nt2]
            raw_beta = float(l[5])
            if float(l[6])**2 >= betathresh:
                nsnps+=1
                upd_pval_beta = float(l[6])
            else:
                upd_pval_beta = 0.0
            rs_id_map[rs_id] = {'chrom':chrom, 'pos':pos, 'nts':nts, 'raw_beta':raw_beta,
                                'upd_pval_beta':upd_pval_beta}
    print "Number of SNPs in block: There are %d SNPs with beta larger than %0.10f"%(nsnps,betathresh)

    return rs_id_map





def parse_ldpred_res_beta2_bounds(file_name,betalow=None,betaup=None):
    rs_id_map = {}
    """
    chrom    pos    sid    nt1    nt2    raw_beta    ldpred_inf_beta    ldpred_beta
    1, 798959, rs11240777, C, T, -1.1901e-02, 3.2443e-03, 2.6821e-04
    """
    nsnps=0
    with open(file_name,'r') as f:
        f.next()
        for line in f:
            l = line.split()
            chrom_str = l[0]
            chrom = int(chrom_str[6:])
            pos = int(l[1])
            rs_id = l[2].strip()
            nt1 = l[3].strip()
            nt2 = l[4].strip()
            nts = [nt1,nt2]
            raw_beta = float(l[5])
            if  float(l[6])**2 >= betalow and float(l[6])**2 < betaup:
                upd_pval_beta = float(l[6])
                nsnps+=1
            else:
                upd_pval_beta = 0.0
            rs_id_map[rs_id] = {'chrom':chrom, 'pos':pos, 'nts':nts, 'raw_beta':raw_beta,
                                'upd_pval_beta':upd_pval_beta}
    print "Number of SNPs in block: There are %d SNPs with beta less than %0.10f and larger than %0.10f  "%(nsnps,betaup, betalow)
    return rs_id_map

def parse_ldpred_res_beta2_low(file_name,betathresh=None):
    rs_id_map = {}
    """
    chrom    pos    sid    nt1    nt2    raw_beta    ldpred_inf_beta    ldpred_beta
    1, 798959, rs11240777, C, T, -1.1901e-02, 3.2443e-03, 2.6821e-04
    """
    nsnps=0
    with open(file_name,'r') as f:
        f.next()
        for line in f:
            l = line.split()
            chrom_str = l[0]
            chrom = int(chrom_str[6:])
            pos = int(l[1])
            rs_id = l[2].strip()
            nt1 = l[3].strip()
            nt2 = l[4].strip()
            nts = [nt1,nt2]
            raw_beta = float(l[5])
            if float(l[6])**2 < betathresh:
                upd_pval_beta = float(l[6])
                nsnps+=1
            else:
                upd_pval_beta = 0.0
            rs_id_map[rs_id] = {'chrom':chrom, 'pos':pos, 'nts':nts, 'raw_beta':raw_beta,
                                'upd_pval_beta':upd_pval_beta}
    print "Number of SNPs in block: There are %d SNPs with beta less than %0.10f"%(nsnps,betathresh)
    return rs_id_map
def parse_ldpred_res(file_name):
    rs_id_map = {}
    """
    chrom    pos    sid    nt1    nt2    raw_beta    ldpred_inf_beta    ldpred_beta
    1, 798959, rs11240777, C, T, -1.1901e-02, 3.2443e-03, 2.6821e-04
    """
    with open(file_name,'r') as f:
        f.next()
        for line in f:
            l = line.split()
            chrom_str = l[0]
            chrom = int(chrom_str[6:])
            pos = int(l[1])
            rs_id = l[2].strip()
            nt1 = l[3].strip()
            nt2 = l[4].strip()
            nts = [nt1,nt2]
            raw_beta = float(l[5])
            upd_pval_beta = float(l[6])
            rs_id_map[rs_id] = {'chrom':chrom, 'pos':pos, 'nts':nts, 'raw_beta':raw_beta,
                                'upd_pval_beta':upd_pval_beta}
    return rs_id_map

def parse_pt_res_lasso(file_name):
    rs_id_map = {}
    """
    chrom    pos    sid    nt1    nt2    raw_beta    raw_pval_beta    upd_beta    upd_pval_beta upd_pval_beta_lasso
    1    798959    rs11240777    C    T    -1.1901e-02    -1.1901e-02    2.6821e-04    2.6821e-04
    """
    with open(file_name,'r') as f:
        f.next()
        for line in f:
            l = line.split()
            chrom_str = l[0]
            chrom = int(chrom_str[6:])
            pos = int(l[1])
            rs_id = l[2].strip()
            nt1 = l[3].strip()
            nt2 = l[4].strip()
            nts = [nt1,nt2]
            raw_beta = float(l[5])
            upd_pval_beta = float(l[8])
            upd_pval_beta_lasso = float(l[9])
            if raw_beta!=0:
                rs_id_map[rs_id] = {'chrom':chrom, 'pos':pos, 'nts':nts, 'raw_beta':raw_beta,
                                    'upd_pval_beta':upd_pval_beta,'upd_pval_beta_lasso':upd_pval_beta_lasso}
                non_zero_chromosomes.add(chrom)

    return rs_id_map

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]



def calc_risk_scores(bimfile_name, rs_id_map, phen_map, K_bins=1,out_file=None,verbose=False):
    num_individs = len(phen_map)
    assert num_individs > 0, 'No individuals found.  Problems parsing the phenotype file?'
    #print K_bins
    if K_bins>1:
        prs_dict_bins={}
        bk = 1
        while bk <= K_bins:
            prs_dict_bins["pval_derived_effects_prs_bin_%d"%bk]=sp.zeros(num_individs)
            bk+=1

    #print prs_dict_bins.keys()
    if bimfile_name is not None:
        raw_effects_prs = sp.zeros(num_individs)
        pval_derived_effects_prs = sp.zeros(num_individs)

        bimf1 = re.sub(r"\[1:22\]", "[0-9]", bimfile_name)
        bimf2 = re.sub(r"\[1:22\]", "[0-2][0-9]", bimfile_name)
        bimfile_list=glob.glob(bimf1+".bim")+glob.glob(bimf2+".bim")
        bimfile_list.sort(key=natural_keys)

        for bimfile in bimfile_list:
            genotype_file = re.sub(r".bim", "", bimfile)
            print 'Get PRS on file %s' % bimfile
            prs_dict = get_prs_bins(genotype_file, rs_id_map, phen_map=phen_map,K_bins=K_bins,verbose=verbose)

            raw_effects_prs += prs_dict['raw_effects_prs']
            pval_derived_effects_prs += prs_dict['pval_derived_effects_prs']
            if K_bins > 1:
                bk = 1
                while bk <= K_bins:
                    prs_dict_bins["pval_derived_effects_prs_bin_%d" % bk] += prs_dict["pval_derived_effects_prs_bin_%d" % bk]
                    bk+=1

        true_phens = prs_dict['true_phens']

    raw_eff_corr = sp.corrcoef(raw_effects_prs, prs_dict['true_phens'])[0, 1]
    raw_eff_r2 = raw_eff_corr ** 2
    pval_eff_corr = sp.corrcoef(pval_derived_effects_prs, prs_dict['true_phens'])[0, 1]
    pval_eff_r2 = pval_eff_corr ** 2

    print 'Final raw effects PRS correlation: %0.4f' % raw_eff_corr
    print 'Final raw effects PRS r2: %0.4f' % raw_eff_r2
    print 'Final LDpred-funct-inf PRS correlation: %0.4f' % pval_eff_corr
    print 'Final LDpred-funct-inf  PRS r2: %0.4f' % pval_eff_r2

    if K_bins > 1:
        X=sp.ones((num_individs, 1))
        bk = 1
        while bk <= K_bins:
            prs_dict_bins["pval_derived_effects_prs_bin_%d" % bk].shape = (len(prs_dict_bins["pval_derived_effects_prs_bin_%d" % bk]), 1)
            X=sp.hstack([prs_dict_bins["pval_derived_effects_prs_bin_%d" % bk], X])
            bk+=1

        true_phens = sp.array(true_phens)
        true_phens.shape = (len(true_phens), 1)

        (betas, rss0, r, s) = linalg.lstsq(X, true_phens)

        ### In sample fit
        Y_pred=sp.dot(X,betas)
        Y_pred.shape=(len(true_phens), )
        # Report prediction accuracy
        bin_in_sample_eff_corr = sp.corrcoef(Y_pred, prs_dict['true_phens'])[0, 1]
        bin_eff_r2 = bin_in_sample_eff_corr ** 2

        print 'Final in-sample LDpredfunct (%d bins) PRS correlation: %0.4f' % (K_bins, bin_in_sample_eff_corr)
        print 'Final in-sample LDpredfunct (%d bins) PRS R2: %0.4f' % (K_bins, bin_eff_r2)
        print 'Final in-sample LDpredfunct (%d bins) PRS adjusted-R2: %0.4f' % (K_bins, 1-(1-bin_eff_r2)*(len(true_phens)-1)/(len(true_phens)-K_bins-1) )

        ###

        test_size= len(true_phens)
        cv_fold_size = int(test_size / 10)
        bound_cv_test = []
        for k in range(10):
            bound_cv_test.append(k * cv_fold_size)
        bound_cv_test.append(test_size - 1)
        bin_eff_r2_arr=[]
        for cv_iter in range(10):
            Xtrain=sp.copy(X)
            Xtest=sp.copy(X)
            Ytrain=sp.copy(true_phens)
            Ytest=sp.copy(true_phens)

            Xtest=Xtest[bound_cv_test[cv_iter]:bound_cv_test[cv_iter+1],]
            Ytest=Ytest[bound_cv_test[cv_iter]:bound_cv_test[cv_iter+1]]

            Xtrain=sp.delete(Xtrain,range(bound_cv_test[cv_iter],bound_cv_test[cv_iter+1]),0)
            Ytrain = sp.delete(Ytrain,range(bound_cv_test[cv_iter],bound_cv_test[cv_iter+1]),0)
            (betas, rss0, r, s) = linalg.lstsq(Xtrain, Ytrain)
            Y_pred = sp.dot(Xtest, betas)
            Y_pred.shape = (len(Ytest),)
            Ytest.shape = (len(Ytest),)
            # Report prediction accuracy
            bin_in_sample_eff_corr = sp.corrcoef(Y_pred, Ytest)[0, 1]
            bin_eff_r2 = bin_in_sample_eff_corr ** 2
            bin_eff_r2_arr.append(bin_eff_r2)

        print 'Final 10-fold cross validation LDpredfunct (%d bins) PRS average R2 : %0.4f ' % (K_bins, sp.mean(bin_eff_r2_arr))



    res_dict = {'pred_r2': pval_eff_r2}

    raw_effects_prs.shape = (len(raw_effects_prs), 1)
    pval_derived_effects_prs.shape = (len(pval_derived_effects_prs), 1)
    true_phens = sp.array(true_phens)
    true_phens.shape = (len(true_phens), 1)

    # Store covariate weights, slope, etc.
    weights_dict = {}


    num_individs = len(prs_dict['pval_derived_effects_prs'])

    # Write PRS out to file.


    if out_file != None:
        with open(out_file, 'w') as f:
            out_str = 'IID, true_phens, raw_effects_prs, pval_derived_effects_prs'
            if K_bins>1:
                Kbins_str= ",".join("Bin_%d"%(1+bin_i)  for bin_i in range(K_bins))
                out_str = out_str + ', ' + Kbins_str
            out_str += '\n'
            f.write(out_str)
            for i in range(num_individs):
                out_str = '%s, %0.6e, %0.6e, %0.6e ' % (
                prs_dict['iids'][i], prs_dict['true_phens'][i], raw_effects_prs[i],
                pval_derived_effects_prs[i])
                bins_prs_ind_i=[]
                if K_bins > 1:
                    bk = 1
                    while bk <= K_bins:
                        bins_prs_ind_i.append(str(prs_dict_bins["pval_derived_effects_prs_bin_%d" % bk][i][0]))
                        bk+=1
                    Kbins_str = ', '.join(map(str, bins_prs_ind_i))
                    out_str = out_str + ', ' + Kbins_str
                out_str += '\n'
                f.write(out_str)

        if weights_dict != None:
            oh5f = h5py.File(out_file + '.weights.hdf5', 'w')
            for k1 in weights_dict.keys():
                kg = oh5f.create_group(k1)
                for k2 in weights_dict[k1]:
                    kg.create_dataset(k2, data=sp.array(weights_dict[k1][k2]))
            oh5f.close()
    return res_dict

def main():
    p_dict = parse_parameters()
    print(p_dict)
    non_zero_chromosomes =set()
    print(p_dict)
    #Parse phenotypes
    if p_dict['pf'] is None:
        if p_dict['vgf'] is not None:
            phen_map = parse_phen_file(p_dict['vgf']+'.fam', 'FAM')
        else:
            raise Exception('Validation phenotypes were not found.')
    else:
        phen_map = parse_phen_file(p_dict['pf'], p_dict['pf_format'])
    iids = set(phen_map.keys())
    num_individs = len(phen_map)
    assert num_individs>0, 'No phenotypes were found!'

    res_dict = {}

        #Plot results?
    if p_dict['res_format'] == 'LDPRED-Funct':


        if p_dict['method'] == "ldscore_solve" or p_dict['method'] == "ldscore":
            weights_file = '%s_LDpred-inf-ldscore.txt' % (p_dict['rf'])
        else:
            weights_file = '%s_LDpred-inf.txt' % (p_dict['rf'])
        print(weights_file)
        if os.path.isfile(weights_file):
            K_bins=int(min(100,math.ceil((0.9*num_individs)/300)))
            if K_bins>1:
                out_file = '%s_validation_LDpred-funct_%d_bins.txt' % (p_dict['out'], K_bins)
            else:
                out_file = '%s_validation_LDpred-funct_inf.txt' % (p_dict['out'])
            rsid_map=parse_ldpred_res_bins(weights_file, K_bins=K_bins)
            calc_risk_scores(p_dict['vgf'], rsid_map, phen_map,K_bins=K_bins, out_file=out_file)


    else:
        raise NotImplementedError('Results file format missing or unknown: %s'%p_dict['res_format'])






if __name__ == '__main__':
    main()
