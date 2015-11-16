vargplvm
========

This repository contains both MATLAB and R code for implementing the Bayesian GP-LVM. The MATLAB code is in the subdirectory vargplvm, the R code in vargplvmR.

For a quick description and <b>sample videos / demos</b> check:<br>
http://git.io/A3Uv
    
Bayesian GP-LVM
===============

The model
---------

The Bayesian GP-LVM (Titsias and Lawrence, 2010) is an extension of the traditional GP-LVM where the latent space is approximately marginalised out in a variational fashion (hence the prefix 'vargplvm').

Let us denote $\mathbf{Y}$ as a matrix of observations (here called
outputs) with dimensions $n \times p$, where $n$ rows correspond to
datapoints and $p$ columns to dimensions. In the latent variable model
(LVM) methodology) we assume that these observations come from a
latent (unobserved or input) space $\mathbf{X}$ in $n\times q$,
$q<<p$. GP-LVM assumes that $\mathbf{Y}$ is generated from
$\mathbf{X}$ using a non-linear mapping with a GP prior. Although this
mapping can be integrated out analytically, the latent variables
cannot. Therefore, GP-LVM is traditionally optimised using MAP,
i.e. by minimizing $-\log(p(\mathbf{Y}|\mathbf{X})p(\mathbf{X}))$.

The Bayesian GP-LVM allows for approximately marginalising out
$\mathbf{X}$, so that optimisation is performed by minimizing
$-\log(p(\mathbf{Y})$. This is done by introducing a variational
distributions $q(\mathbf{X}) \sim \mathcal{N}(\mu, \Sigma)$ for which $\mu$ and
$\Sigma$ are to be determined in optimisation (and are variational
parameters). This has some advantages:

a. As a by-product we also obtain the approximate posterior by
  $q(\mathbf{X}) \sim p(\mathbf{X}|\mathbf{Y})$.
b. Optimising using p(\mathbf{Y}) is more robust to overfitting, since we now have
   a proper distribution for the latent variables (rather than point estimates)
c. We can perform automatic dimensionality detection by using Automatic
   Relevance Determination (ARD) covariance functions for the GP prior.
dd ...

In vargplvm, we can think as the means mu to replace the estimates for $\mathbf{X}$.

In the standard Bayesian GP-LVM, $p(\mathbf{X})$ is given a standard
normal prior and is integrated out.

* Variational Gaussian process dynamical systems (VGPDS), Damianou et. al, 2011

What happens when we have some prior information regarding the nature
of our observations $\mathbf{Y}$? For example, if $\mathbf{Y}$ form a
time-series, it is natural to want to constrain $\mathbf{X}$ to also
form a smooth paths. In this case, if $\mathbf{Y}$ is associated with
a time vector t, then we also want $\mathbf{X}$ to be associated
(constrained) by this vector as well.

You can also think as t being an *observed input* to the whole process.
This can be a temporal vector, class labels or whatever. Thus, the VGPDS
extension can be seen as a supervised version of the previous model, the
Bayesian GP-LVM. The code is written with a notation which assumes temporal
t, but you can think of t being any sort of multi-dimensional input without
changing a single line in the code.

Formally, if $\mathbf{Y} \in n \times p$, then $t \in n \times 1$ (or
$n \times p_t$, if we are not talking about time). VGPDS allows the
latent layer $\mathbf{X}$ to have itself a prior which is also a GP
(in contrast to a standard normal prior used in Bayesian GP-LVM). So:
\[\mathbf{Y} = f(\mathbf{X}) + e_1,
\]
where e is Gaussian noise and
\[
f ~ GP(0, K_x(\mathbf{X},\mathbf{X}))
\]
\[
x ~ GP(0, k_t(t,t)) 
\]

and $K_x$ and $K_t$ are constructed from covariance functions $k_x$
and $k_t$.  Although the Baysian GP-LVM framework constrains us to use
a specific set of covariance functions $k_x$ for $\mathbf{X}$ (since
$\mathbf{X}$ is integrated out), we can use any kind of covariance
function for $k_t$. This gives rise to very powerful models, for
example if we know there is some periodicity in the temporal inputs
$t$ we can choose $k_t$ to be a temporal covariance function.

* Reparameterization: In contrast to the $q$ diagonal covariance
matrices Sigma of $q(\mathbf{X}) \sim \mathcal{N}(\mu, \Sigma)$ of Baysian
GP-LVM, in VGPDS each of the $q$ $Sigma_j$ is full-rank $n\times n$
making optimisation very challenging. In VGPDS a reparameterization
trick is used so that each $\Sigma_i$ can be reparameterized with a
vector $\lambda_i$ which is $n \times 1$ (equivalently, a diagonal
matrix $\Lambda_i$).

* Initialising the reparameterized values It turns out (please consult
paper) that $\Lambda_i$ is depends on $\Sigma_i$ in a complex
way. Since we now have access only to $\Lambda_i$ as free parameters
and not to $\Sigma_i$ directly, it is difficult to initialise
$\Sigma_i$ to our desired values (ideally, $\Sigma_i$ would take
initial values in the range of 0.1 to 0.5). The code implements a
simple brute force method which picks a reasonable initial
$\Lambda_i$.


MATLAB Code
===========

Matlab code contains implementations of:
 - (1) the Bayesian GPLVM (see Titsias and Lawrence, AISTATS 2010)
 - (2) Variational Gaussian Process Dynamical Systems (VGPDS) (see Damianou et al., NIPS 2011)
 - (3) Manifold Relevance Determination (MRD) (see Damianou et al., ICML 2012)

Files for (1) and (2) have the prefix "vargplvm". Files for (3) have the prefix "svargplvm".
The rest of the files are generic and used for all methods.

Dependencies graph:
1. GPmat - The GPmat toolbox: https://github.com/SheffieldML/GPmat
2. Netlab v.3.3: http://www1.aston.ac.uk/ncrg/
3. Isomap.m: http://web.mit.edu/cocosci/isomap/code/Isomap.m
4. L2_distance.m: http://web.mit.edu/cocosci/isomap/code/L2_distance.m

L2_distance ---- Isomap ---- GPmat ---- vargplvm
							/			/
				Netlab ----------------

Getting started with the code:
------------------------------
 - Please check vargplvm/html/index.html for a short overview of this package.
 - Check vargplvm/matlab/vargplvm_README.txt for Bayesian GP-LVM and VGPDS.
 - Check vargplvm/matlab/VGPDStutorial.m for introductory demonstrations for VGPDS (see demosDynamics.m for more).
 - Check vargplvm/matlab/svargplvm_README.txt for an overview of MRD
 - Check vargplvm/matlab/MRDtutorial.m for introductory demonstrations for MRD.

R Code
======

The R code was initially implemented by Jie Hao <j.hao@ic.ac.uk> based on the MATLAB code. Further work on the code and documentation was performed by Nicolas Durrande <nicolas.durrande@emse.fr>. Both were funded by the <a href="http://staffwww.dcs.sheffield.ac.uk/people/N.Lawrence/projects/biopredyn/"> EU FP7-KBBE Project Ref 289434 "From Data to Models: New Bioinformatics Methods and Tools for Data-Driven Predictive Dynamic Modelling in Biotechnological Applications"</a>.
