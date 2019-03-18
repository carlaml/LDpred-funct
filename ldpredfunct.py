#!/usr/bin/env python
# Usar ldpred-master/coord_genotypes_ldpredfunct.py como base

#### bim_file = "/n/scratch2/sg374/UKB_LIKELIHOOD/25K_Carla/plink_files/Final_UK10K.[1:22]"

### input plink files: per chromsom
import getopt
import sys
import os
import math
import validate_prs_ldpredfunct as prs
import coord_genotypes_ldpredfunct_v1_2 as coord
import LDpredfunct_bayes_shrink as ldpredfunct
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


    long_options_list = ['FUNCT_FILE=','gf=', 'gmdir=', 'check_mafs', 'coord=', 'maf=', 'skip_coordination','verbose', 'skip_ambiguous',"chisq", 'ssf=', 'N=',"K=", "posterior_means=", 'ld_radius=', 'H2=', 'out=',"pf="]

    p_dict = {'FUNCT_FILE':None,'gf':None, 'gmdir':None, 'check_mafs':False, "coord":"output-coordinated", 'maf':0.01,'K':None,"chisq":False,  'skip_coordination':False, 'skip_ambiguous':False,
    'ssf':None, 'N':None, "posterior_means":"output-posterior_means", 'ld_radius':None, 'H2':None,'verbose':False, 'out':"output-prs","pf":None}

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
            elif opt in ("--gmdir"):
                p_dict['gmdir'] = arg
            elif opt in ("--check_mafs"):
                p_dict['check_mafs'] = True
            elif opt in ("--coord"):
                p_dict['coord'] = arg
            elif opt in ("--skip_coordination"):
                p_dict['skip_coordination'] = True
            elif opt in ("--skip_ambiguous"):
                p_dict['skip_ambiguous'] = True
            elif opt in ("--verbose"):
                p_dict['verbose'] = True
            elif opt in ("--chisq"):
                p_dict['chisq'] = True
            elif opt in ("--ssf"):
                p_dict['ssf'] = arg
            elif opt in ("--pf"):
                p_dict['pf'] = arg
            elif opt in ("--K"):
                p_dict['K'] = arg
            elif opt in ("--FUNCT_FILE"):
                p_dict['FUNCT_FILE'] = arg
            elif opt in ("--N"):
                p_dict['N'] = int(arg)
            elif opt in ("--posterior_means"): p_dict['posterior_means'] = arg
            elif opt in ("--ld_radius"): p_dict['ld_radius'] = int(arg)
            elif opt in ("--H2"): p_dict['H2'] = float(arg)
            elif opt in ("--out"): p_dict['out'] = arg
            else:
                print "Unkown option:", opt
                print "Use -h option for usage information."
                sys.exit(2)
    else:
        print __doc__
        sys.exit(0)
    return p_dict


def main():
    p_dict = parse_parameters()
    print(p_dict)


    if p_dict['N'] is None:
        print 'Please specify an integer value for the sample size used to calculate the GWAS summary statistics.'

    if os.path.isfile(p_dict['out']):
        print 'Output file (%s) already exists!  Delete, rename it, or use a different output file.'%(p_dict['out'])
        raise Exception('Output file already exists!')
    """
        --ssf= Summary statistics file. It most have the following fields (order doesnt matter): CHR SNP BP Z BETA P
        --gf= plink file (one file per-chromosome). Example:/plink_files/Final_UKBiobank.[1:22]
        --coord= Output for hdf5 coordinated genotype and summary statistics file
        --n= Sample size of training
        --FUNCT_FILE= File with heritability enrichemnts. File should have two columns: SNP BaselineLD. Recommendation use the output from S-LDSC under the baselineLD model. See
        gmdir
        check_mafs
        maf
        skip_coordination
        skip_ambiguous
        posterior_means
        ld_radius
        N
        H2
        pf
        out
    """
    print("Step 1: Coordinate summary statistics, genotype and functional enrichments files.\n")
    coord.parse_sum_stats_standard_ldscore(filename=p_dict['ssf'], bimfile_name=p_dict['gf'], hdf5_file_name=p_dict['coord'],
                                     n=p_dict['N'],
                                     outfile=p_dict['coord'] + "_snps_NaN.txt", FUNCT_FILE=p_dict["FUNCT_FILE"],CHISQ=p_dict['chisq'])

    coord.coordinate_genot_ss(genotype_filename=p_dict['gf'], genetic_map_dir=p_dict['gmdir'], check_mafs=p_dict['check_mafs'],
                    hdf5_file_name=p_dict['coord'], min_maf=p_dict['maf'], skip_coordination=p_dict['skip_coordination'],
                    method="STANDARD_FUNCT", skip_ambiguous=p_dict['skip_ambiguous'])

    print("Step 2: Compute posterior mean effect sizes using a functional prior.")
    ldpredfunct.ldpred_funct_genomewide(data_file=p_dict['coord'], out_file_prefix=p_dict['posterior_means'], ld_radius=p_dict['ld_radius'],
                          n=p_dict['N'],h2=p_dict['H2'],verbose=p_dict['verbose'])

    print("Step 3: Compute polygenic risk score using previously computed posterior mean effect sizes.")

    if p_dict['pf'] is None:
        if p_dict['gf'] is not None:
            phen_map = prs.parse_phen_file(p_dict['gf']+'.fam', 'FAM')
        else:
            raise Exception('Validation phenotypes were not found.')
    else:
        phen_map = prs.parse_phen_file(p_dict['pf'],'STANDARD')
    iids = set(phen_map.keys())
    num_individs = len(phen_map)
    assert num_individs>0, 'No phenotypes were found!'

    weights_file = '%s_LDpred-inf-ldscore.txt' % (p_dict['posterior_means'])
    print("Reading %s"%weights_file)
    if os.path.isfile(weights_file):
        if p_dict['K']==None:
            K_bins=int(min(100,math.ceil((0.9*num_individs)/300)))
        else:
            K_bins=int(p_dict['K'])
        if K_bins>1:
            out_file = '%s_validation_LDpred-funct_%d_bins.txt' % (p_dict['out'], K_bins)
        else:
            out_file = '%s_validation_LDpred-funct_inf.txt' % (p_dict['out'])
        rsid_map=prs.parse_ldpred_res_bins(weights_file, K_bins=K_bins)
        prs.calc_risk_scores(p_dict['gf'], rsid_map, phen_map,K_bins=K_bins, out_file=out_file,verbose=p_dict['verbose'])


if __name__ == '__main__':
    main()

