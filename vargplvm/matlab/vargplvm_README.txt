% VARGPLVM_README
% Copyright: Andreas C. Damianou, 2012, 2013
% VARGPLVM
                             _                 
                            | |                
 __   ____ _ _ __ __ _ _ __ | |_   ___ __ ___  
 \ \ / / _` | '__/ _` | '_ \| \ \ / / '_ ` _ \ 
  \ V / (_| | | | (_| | |_) | |\ V /| | | | | |
   \_/ \__,_|_|  \__, | .__/|_| \_/ |_| |_| |_|
                  __/ | |                      
                 |___/|_|                      




                              R E A D M E

___________________________________________________________________________
##################### 1. DESCRIPTION AND CHANGELOG ########################
___________________________________________________________________________


The vargplvm toolbox is an implementation of the variational approximation for the Bayesian GPLVM.

Version 0.3
-----------

The software is now merged with MRD (svargplvm).

!! NOTE: If you are interested in MRD, please check svargplvm_README.txt.


Version 0.2
-----------

The software now is extended with dynamics.


Version 0.1
-----------

First version of the software with implementation of the 2010 AISTATS paper.


___________________________________________________________________________
##################### 2. DEPENDENCIES GRAPH ########################
___________________________________________________________________________

(1) GPmat - Neil Lawrence's GP matlab toolbox: https://github.com/SheffieldML/GPmat
(2) Netlab v.3.3: http://www1.aston.ac.uk/ncrg/
(3) Isomap.m: http://web.mit.edu/cocosci/isomap/code/Isomap.m
(4) L2_distance.m: http://web.mit.edu/cocosci/isomap/code/L2_distance.m

L2_distance ---- Isomap ---- GPmat ---- vargplvm
							/			/
				Netlab ----------------


___________________________________________________________________________
##################### 3. GENERAL ########################
___________________________________________________________________________

_____________________________ The model (theory)

# Bayesian GP-LVM (Titsias and Lawrence, 2010):

The Bayesian GP-LVM is an extension of the traditional GP-LVM where the
latent space is approximately marginalised out in a variational fashion
(hence the prefix 'vargplvm').

Let us denote Y as a matrix of observations (here called outputs) with
dimensions N x D, where N rows correspond to datapoints and D columns to
dimensions. In the latent variable model (LVM) methodology) we assume that
these observations come from a latent (unobserved or input) space X in
N x Q, Q << D. GP-LVM assumes that Y is generated from X using a non-linear
mapping with a GP prior. Although this mapping can be integrated out
analytically, the latent variables cannot. Therefore, GP-LVM is traditionally
optimised using MAP, i.e. by minimizing -log(p(Y|X)p(X)).

The Bayesian GP-LVM allows for approximately marginalising out X, so that
optimisation is performed by minimizing -log(p(Y). This is done by introducing
a variational distributions q(X) ~ N(mu, Sigma) for which mu and Sigma are
to be determined in optimisation (and are variational parameters). This has
some advantages:
a) As a by-product we also obtain the approximate posterior by
  q(X) ~ p(X|Y).
b) Optimising using p(Y) is more robust to overfitting, since we now have
   a proper distribution for the latent variables (rather than point estimates)
c) We can perform automatic dimensionality detection by using Automatic
   Relevance Determination (ARD) covariance functions for the GP prior.
d) ...

In vargplvm, we can think as the means mu to replace the estimates for X.

In the standard Bayesian GP-LVM, p(X) is given a standard normal prior
and is integrated out.

# Variational Gaussian process dynamical systems (VGPDS), Damianou et. al, 2011

What happens when we have some prior information regarding the nature of
our observations Y? For example, if Y form a time-series, it is natural to
want to constrain X to also form a smooth paths. In this case, if Y is
associated with a time vector t, then we also want X to be associated
(constrained) by this vector as well. 

You can also think as t being an *observed input* to the whole process.
This can be a temporal vector, class labels or whatever. Thus, the VGPDS
extension can be seen as a supervised version of the previous model, the
Bayesian GP-LVM. The code is written with a notation which assumes temporal
t, but you can think of t being any sort of multi-dimensional input without
changing a single line in the code.

Formally, if Y \in N x D, then t \in N x 1 (or N x D_t, if we are not talking
about time). VGPDS allows the latent layer X to have itself a prior which
is also a GP (in contrast to a standard normal prior used in Bayesian
GP-LVM). So:
Y = f(X) + e_1, where e is Gaussian noise and
f ~ GP(0, K_x(X,X))
x ~ GP(0, K_t(t,t)) 

and K_x and K_t are constructed from covariance functions k_x and k_t.
Although the Baysian GP-LVM framework constrains us to use a specific set
of covariance functions k_x for X (since X is integrated out), we can use
any kind of covariance function for k_t. This gives rise to very powerful
models, for example if we know there is some periodicity in the temporal
inputs t we can choose k_t to be a temporal cov. function.

* Reparameterization:
In contrast to the Q diagonal covariance matrices Sigma of q(X) ~ N(mu, Sigma)
of Baysian GP-LVM, in VGPDS each of the Q Sigma_q is full-rank N x N and makes optimisation very
challenging. In VGPDS a reparameterization trick is used so that each Sigma_q
can be reparameterized with a vector lambda_q which is N x 1 (equivalently,
a diagonal matrix Lambda_q).

* Initialising the reparameterized values
It turns out (please consult paper) that Lambda_q is depending on Sigma_q in
a complex way. Since we now have access only to Lambda_q as free parameters
and not to Sigma_q directly, it is difficult to initialise Sigma_q to
our desired values (ideally, Sigma_q would take initial values in the range
of 0.1 to 0.5). The code implements a simple brute force method which picks
a reasonable initial Lambda_q.



_____________________________ The model (in the code)

...

_____________________________ Demos and global configuration to run experiments
...

_____________________________ Initialisation

Check the vargplvmModelCreate in combination with vargplvm_init.
Also check "* Initialising the reparameterized values" section above.

_____________________________ Optimisation

The optimisation can be done while initialising the variational
distribution by keeping the SNR fixed (by controlling model.beta to have a
constant ratio w.r.t the variance of the data -which is equal to
model.kern.sigmaf that is also fixed, when the 'rbfardjit' kernel is
used).

In the high level demos (i.e. the ones that optimise the model using the
function vargplvmOptimiseModel.m rather than directly the vargplvmOptimise.m),
it is enough to just set the desired iterations for initialising the var.
distr (initVardistIters = ...) and the initial SNR (initSNR, as above).
Otherwise this has to be done manually by doing one optimisation round with
keeping beta fixed (control model.learnBeta and model.learnSigmaf) and one
optimisation round by releasing it. Check the demos, in both cases there are
examples.

In any case, it's a good practise to call vargplvmShowSNR(model) to make sure
that the SNR of the optimised model is relatively high (e.g. at least 10).

!! Note that if rbfardjit is used as a mapping kernel, then its variance
(sigmaf) is always initialised to the variance of the data, var(m(:)),
so by keeping both sigmaf and beta fixed the SNR is fixed in a better
way that just fixing beta and using some other mapping kernel. 


___________________ Avoiding bad local optima

The optimisation procedure is gradient based, ie there is no an analytic
form to a unique solution. This means that the quality of the optimisation
is depending upon many things, mostly: initialisation, optimiser used,
numerical errors.

From the above, the easiest to control is the initialisation.
After optimisation, call svargplvmCheckSNR(vargplvmShowSNR(model))
to check if a bad local minimum is reached (the optimisation function actually
calls automatically this check).
If this a bad local optimum is reached (ie signal is very low compared to
noise), then the model is badly trained and the results are unreliable.

When this unfortunate scenario happens, please try a different initialisation,
for example:
    - increase the number of iterations for init. the var. distr. (initVardistIters).
    - increase the initial SNR
    - set the reparametrized parameters Lambda for Sigma (code notation: model.dynamics.vardist.covars)
      to a value that results in more reasonable initial Sigma (i.e. pick a Lambda_q which results
      in more reasonable Sigma_q where more reasonable usually means smaller). 
      See the above section: "* Initialising the reparameterized values".
    - initialising the means of q(X) with a different method than the default 'pca',
      e.g. with the outputs Y themselves (if D == Q) or a 'gplvm', etc.
    - preprocess the data first, so that some noise is removed.
    - ...


___________________________________________________________________________
##################### 4. References ########################
___________________________________________________________________________

M. Titsias, N. Lawrence, Bayesian Gaussian process latent variable model, 
AISTATS 2010

A. Damianou, M. Titsias, N. Lawrence, Variational Gaussian process
dynamical systems, NIPS 2011

___________________________________________________________________________
################ 5. Manifold Relevance Determination (MRD) ################
___________________________________________________________________________

Please consult svargplvm_README.txt
