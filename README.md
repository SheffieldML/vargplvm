vargplvm
========

Thsi repository contains both MATLAB and R code for implementing the Bayesian GP-LVM. The MATLAB code is in the subdirectory vargplvm, the R code in vargplvmR.

MATLAB Code
===========

Matlab code contains implementations of:
 - (1) the Bayesian GPLVM (see Titsias and Lawrence, AISTATS 2010)
 - (2) Variational Gaussian Process Dynamical Systems (VGPDS) (see Damianou et al., NIPS 2011)
 - (3) Manifold Relevance Determination (MRD) (see Damianou et al., ICML 2012)

Files for (1) and (2) have the prefix "vargplvm". Files for (3) have the prefix "svargplvm".
The rest of the files are generic and used for all methods.

Dependencies graph:
(1) GPmat - The GPmat toolbox: https://github.com/SheffieldML/GPmat
(2) Netlab v.3.3: http://www1.aston.ac.uk/ncrg/
(3) Isomap.m: http://web.mit.edu/cocosci/isomap/code/Isomap.m
(4) L2_distance.m: http://web.mit.edu/cocosci/isomap/code/L2_distance.m

L2_distance ---- Isomap ---- GPmat ---- vargplvm
							/			/
				Netlab ----------------

Getting started:
 - Please check vargplvm/html/index.html for a short overview of this package.
 - Check vargplvm/matlab/vargplvm_README.txt for Bayesian GP-LVM and VGPDS.
 - Check vargplvm/matlab/demosDynamics.m for introductory demonstrations for VGPDS.
 - Check vargplvm/matlab/svargplvm_README.txt for an overview of MRD
 - Check vargplvm/matlab/MRDtutorial.m for introductory demonstrations for MRD.

R Code
======

The R code was initially implemented by Jie Hao <j.hao@ic.ac.uk> based on the MATLAB code. Further work on the code and documentation was performed by Nicolas Durrande <nicolas.durrande@emse.fr>. Both were funded by the <a href="http://staffwww.dcs.sheffield.ac.uk/people/N.Lawrence/projects/biopredyn/"> EU FP7-KBBE Project Ref 289434 "From Data to Models: New Bioinformatics Methods and Tools for Data-Driven Predictive Dynamic Modelling in Biotechnological Applications"</a>.