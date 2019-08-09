# Author: Carla Marquez-Luna 
# For full documentation: https://github.com/carlaml/LDpred-funct
# You need to install the packages h5py, scipy, and libplinkio. With pip (see https://pip.pypa.io/en/latest/quickstart.html), one can install libplinkio using the following command: pip install plinkio
# Running time in a PC for test data: < 2 minutes 

# Input files:
plinkfile="test/TSI_[1:22]"
phenotype="test/TSI_simulated_trait.txt"
statsfile="test/summary_statistics_traininig.txt"
functfile=test/test_functfile.txt

# Output files
outCoord="test/Coord_Final"
outLdpredfunct="test/ldpredfunct_posterior_means"
outValidate="test/ldpredfunct_prs"

N=100
h2=1

# Note: Make sure that the file ${outCoord}  does not exists already. 
# rm ${outCoord} 
python LDpred-funct-master/ldpredfunct.py --gf=${plinkfile} --pf=${phenotype} --FUNCT_FILE=${functfile}  --coord=${outCoord} --ssf=${statsfile} --N=${N} --posterior_means=${outLdpredfunct}  --H2=${h2} --out=${outValidate} > ${outValidate}.log
