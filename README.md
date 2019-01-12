# LDpred-funct
Here we present the software for the method LDpredfunct described in Marquez-Luna, et al, "Modeling functional enrichment improves polygenic prediction accuracy in UK Biobank and 23andMe data sets", Biorxiv. https://www.biorxiv.org/content/early/2018/07/24/375337 

### Installing LDpred-funct

As with most Python packages, configurating LDpred-funct is simple. You can either use git (which is installed on most systems) and clone this repository using the following git command:

git clone https://github.com/carlaml/LDpred-funct.git

Alternatively, you can simply download the source files and place them somewhere.

You also need to install the packages **h5py**, **scipy**, and **libplinkio**. With pip (see https://pip.pypa.io/en/latest/quickstart.html), one can install **libplinkio** using the following command: ```pip install plinkio```


### Input files

1. Plink files. **One plink file (binary PED format) per chromosome** from the validation (insert the character "[1:22]" instead of the chromsome numbers). These are used as LD-reference panel and later to compute PRS.

  - Example: plinkfile="my_plink_files/Final.chrom[1:22]"
  - Use flag: --gf
 
2. File with functional enrichments (see Gazal et al 2017 Nat Genet and Finucane et al 2015 Nat Genet).
    - First you will need to estimate the per-SNP heritability inferred using S-LDSC (see instructions [here](https://github.com/bulik/ldsc/wiki/Partitioned-Heritability)) under the baselineLD model (you can download Baseline-LD annotations [here](https://data.broadinstitute.org/alkesgroup/LDSCORE/)). When running S-LDSC make sure to use the --print-coefficients flag to get the regression coefficients.
    - After running S-LDSC:
      1. Get h2g estimate from the *.log file.
      2. Get the regression coefficients from the *.results file (column 8). Divide the regression coeffients by h2g, define it as T=tau/h2g which is a vector of dimension Cx1, where C is the total number of annotations.
      3. From the baselineLD annotations downloaded from [here](https://data.broadinstitute.org/alkesgroup/LDSCORE/), read the annotations file baselineLD.*.annot.gz, and only keep the annotations columns (i.e. remove first 4 columns). Call this matrix X, with dimensions MxC, where M is the number of SNPs and C is the total number of annotations.
      4. Define the expected per-SNP heritability (sigma2_i from the LDpred-funct manuscript) as the result from multiplying the X and T.

  - Format of FUNCTFILE:
    - Column 1: SNP ID
    - Column 2: per-SNP heritability

  - Use flag: --FUNCT_FILE

3. Summary statistics file. Please check that the summary statistics contains a column for each of the following field (header is important here, important fields are highlighted in bold font, the order of the columns it is not important).
    - **CHR**   Chromosome
    - **SNP**   SNP ID
    - **BP**    Physical position (base-pair)
    - **A1**    Minor allele name (based on whole sample) 
    - **A2**    Major allele name 
    - **P**     Asymptotic p-value  
    - **BETA**  Effect size
    - **Z**     Z-score (default). If instead of Z-score the Chi-square statistic is provided, use the flag --chisq, and **CHISQ** as column field.
    

  - Use flag: --ssf

4. Phenotype file.
  - Format:
    - Column 1: FID
    - Column 2: phenotype
  - This file doesn't have a header.

  - Use flag: --pf


### Input parameters:
1. Training sample size. 
  - Use flag: --N
2. Estimated SNP Heritability (pre-compute this using your favorite method). 
  - Use flag: --H2
3. LD radius (optional). If not provided, it is computed as (1/2)*0.15% of total number of SNPs.  
  - Use flag: --ld_radius

### Output files
1. Coordinated files: This is an hdf5 file that stores the coordinated genotype data with the summary statistics and functional enrichments.
  - Use flag: --coord
  - Note: the output file needs to be named differently for different runs.
  
2. Posterior mean effect sizes: Estimated posterior mean effect size from LDpred-funct-inf.
  - Use flag: --posterior_means
  
3. Output: Polygenic risk score for each individual in the validation. 
  - Description:
    - Column 1: Sample ID
    - Column 2: True phenotype
    - Column 3: PRS using all-snps and marginal effect sizes.
    - Colunm 4: PRS obtained using LD-pred-funct-inf
    - Column 5-K: PRS(k) defined in equation 5 from Marquez-Luna, et al, Biorxiv.

  - Use flag: --out
  
### Example
```
plinkfile="my_plink_files/Final.chrom[1:22]"
outCoord="my_analysis/Coord_Final"
statsfile="my_analysis/sumary_statistics.txt"
outLdpredfunct="my_analysis/ldpredfunct_posterior_means"
outValidate="my_analysis/ldpredfunct_prs"
phenotype="my_analysis/trait_1.txt"
N=training sample size
h2= pre-computed SNP heritability

optional flags:
--ld_radius= Pre-defined ld-radius
--K= Number of bins for LDpred-funct

python LDpred-funct/ldpredfunct.py --gf=${plinkfile} --pf=${phenotype} --FUNCT_FILE=${statsfile}  --coord=${outCoord} --ssf=${statsfile} --N=${N} --posterior_means=${outLdpredfunct}  --H2=${h2} --out=${outValidate} > ${outValidate}.log
```
### Running recommendations
For UK Biobank analysis described in Marquez-Luna, et al, "Modeling functional enrichment improves polygenic prediction accuracy in UK Biobank and 23andMe data sets", Biorxiv. https://www.biorxiv.org/content/early/2018/07/24/375337, I requested 70G of memory and the jobs lasted approximately 10 hours.  

#### Acknowledgments
We used LDpred software developed by Bjarni Vilhjalmsson as basis for this software.
