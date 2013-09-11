vargplvm
========

Matlab code for:
 - (1) the Bayesian GPLVM (see Titsias and Lawrence, AISTATS 2010)
 - (2) Variational Gaussian Process Dynamical Systems (VGPDS) (see Damianou et al., NIPS 2011)
 - (3) Manifold Relevance Determination (MRD) (see Damianou et al., ICML 2012)

Files for (1) and (2) have the prefix "vargplvm". Files for (3) have the prefix "svargplvm".
Rest of the files are generic and used for all methods.

Dependencies graph:
- (1) GPmat - Neil Lawrence's GP matlab toolbox: https://github.com/SheffieldML/GPmat
- (2) Netlab v.3.3: http://www1.aston.ac.uk/ncrg/
- (3) Isomap.m: http://web.mit.edu/cocosci/isomap/code/Isomap.m
- (4) L2_distance.m: http://web.mit.edu/cocosci/isomap/code/L2_distance.m
- (5) keep.m: http://www.mathworks.com/matlabcentral/fileexchange/181-keep/content/keep.m

L2_distance ---- Isomap ---- GPmat ---- vargplvm
							/			/
				Netlab ----------------

Getting started:
 - Please check vargplvm/html/index.html for a short overview of this package.
 - Check vargplvm/matlab/vargplvm_README.txt for Bayesian GP-LVM and VGPDS.
 - Check vargplv/matlab/demosDynamics.m for introductory demonstrations for VGPDS.
 - Check vargplvm/matlab/svargplvm_README.txt for an overview of MRD
 - Check vargplvm/matlab/MRDtutorial.m for introductory demonstrations for MRD.
