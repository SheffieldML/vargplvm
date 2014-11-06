% VARGPLVM_README
% Copyright: Andreas C. Damianou, 2012, 2013, 2014
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

Version 1.0
-----------

The current release includes the Manifold Relevance Determination (MRD) 
method (shared variational GP-LVM) and various utilities for easier
initialisation and demo configuration, as well as better documentation.



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










_____________________________ The model (in the code)_____________________



% THE MODEL STRUCTURE FOR the vargplvm.
% Copyright: Andreas C. Damianou, 2014
%

%------ VARGPLVM MODEL (basic structure) ----------%

%- Example of an instantiation from the tutorial
>> model
model = 

         type: 'vargplvm'     % The name of the model
       approx: 'dtcvar'
  learnScales: 0
 optimiseBeta: 1
betaTransform: 'exp'
            q: 3
            d: 5
            N: 50
    optimiser: 'scg2'
         bias: [0.0645 0.6978 0.6973 -0.3585 0.2118]
        scale: [1 1 1 1 1]
            y: [50x5 double]
            m: [50x5 double]
     KLweight: 0.5000
            X: [50x3 double]
    learnBeta: 1
         DgtN: 0
         TrYY: 175.2405
         date: '10-Oct-2014'
         kern: [1x1 struct]
      vardist: [1x1 struct]
            k: 20
  fixInducing: 0
          X_u: [20x3 double]
         beta: 142.0905
        prior: [1x1 struct]
  initVardist: 1
  learnSigmaf: 0
    numParams: 365
      nParams: 369
         K_uu: [20x20 double]
         Psi0: 35.1889
         Psi1: [50x20 double]
         Psi2: [20x20 double]
           Lm: [20x20 double]
        invLm: [20x20 double]
       invLmT: [20x20 double]
            C: [20x20 double]
          TrC: 35.1811
           At: [20x20 double]
          Lat: [20x20 double]
       invLat: [20x20 double]
      invLatT: [20x20 double]
     logDetAt: -52.2889
           P1: [20x20 double]
            P: [20x5 double]
         TrPP: 173.7928
            B: [20x5 double]
      invK_uu: [20x20 double]
           T1: [20x20 double]
  dataSetInfo: [1x1 struct]
    globalOpt: [1x1 struct]
     dynamics: [1x1 struct]


%-- Explanation of the fields

         type: The name of the model. This description is for 'vargplvm'.
       approx: Inducing point approximation used. This model only uses 'dtcvar', so you can't change it.
               'dtcvar' comes from Titsias 09 and Titsias and Lawrence 2010 paper.
  learnScales: Don't worry about that field (leave default value).
 optimiseBeta: Don't worry about that field (leave default value).
betaTransform: Don't worry about that field (leave default value).
            q: The dimensionality of X.
            d: The dimensionality of Y.
            N: The number of training points.
    optimiser: Which optimiser to use for training the model (i.e. for minimising the objective function
               given the gradients function --check vargplvmOptimise.m). 'scg2' is the default.
         bias: GPs have to be zero mean. Therefore, your observed data Y must first be centered.
               In vargplvmCreate.m the mean(model.y) is computed and stored in model.bias. This value is
               then subtracted by the data to center them, but must be stored so that in predictions it's
               added back to get to the original space.
        scale: Apart from centering the data you might also want to scale them. This is not required
               but if you select to do it (via the options passed to vargplvmCreate.m) then again you
               need to store the value that you divide with, so that you can multiply back in predictions
               (or when you want your original data).
            y: This is the raw data you pass to vargplvmCreate.m, it's the observed outputs.
            m: This is the data that is actually used by the model internally, it's model.y after bias
               and scale (see above) have been applied.
     KLweight: Don't worry about that field and you're advised to leave it to the default value 0.5.
            X: Your latent space. It is a copied version of model.vardist.means.
    learnBeta: Depricated (see initVardist).
         DgtN: If model.d >> model.n, then there's a trick you can do to speed up computations.
               Read "Extensions" in Damianou et al. 2011, NIPS. This field here says if this "trick"
               is activated. If it is, then model.m is replaced by cholesky(model.m * model.m').
               Since the data only appear in the form tr(Y Y^T), the above substitution does not
               affect computations but makes them much faster, since model.m now is much smaller in size.
         TrYY: trace of model.m*model.m'
         date: Date of the model created
         kern: The kernel of the mapping between X and Y. See below.
      vardist: The variatinal distribution structure representing q(X). See below.
            k: The number of inducing points.
  fixInducing: If true, then you need k == model.n, and the inducing points X_u will be tied
               to model.X (ie they're going to be represented as 1 set of parameters and their
               gradients will be added). The intuition behind this is that, the optimal solution
               is to exactly have X_u == X, but usually we have k << n.
          X_u: The matrix of inducing points.
         beta: The noise variance, ie y = f(x) + e, e~N(0, beta^{-1}I).
        prior: Don't worry about that field (leave default value).
  initVardist: When optimising the model, it's wise to first fix some of the very fluctuating parameters
               and optimise the rest by keeping the aforementioned set fixed to their initialisations.
               Specifically, the fixed parameters are beta and the variance of the rbfardjit kernel, 
               which correspond to the Signal to Noise Ratio which initially has to be forced to be high
               to avoid bad local minima close to bad initialisations.
  learnSigmaf: Depricated (see initVardist).
    numParams: The number of parameters of the model.
      nParams: Depricated but unfortunately still used (refactoring not done yet and there's some 
               inconsistency that fortunatelly doesn't affect the code due to workarounds).
         K_uu: The covariance matrix evaluated on the inducing inputs X_u.
         Psi0: Psi0 statistic (see BGPLVM paper)
         Psi1: Psi1 statistic (see BGPLVM paper)
         Psi2: Psi2 statistic (see BGPLVM paper)
  %--- The following parameters are just quantities used by the model to evaluate the objective and
       the gradients. Because the objective and gradients are in different functions, we dont want to
       compute common quantities twice, so we store them in the model structure. Generally you shouldn't
       worry about the values of the following fields, but it might be useful to explore them when you
       debug (e.g. At becomes badly conditioned, or e.g. C becomes huge and causes numerical
       instabilities)
  %------------------------------
           Lm: 
        invLm: 
       invLmT: 
            C: 
          TrC: 
           At: 
          Lat: 
       invLat: 
      invLatT: 
     logDetAt: 
           P1: 
            P: 
         TrPP: 
            B: 
      invK_uu: 
           T1: 
    %---- End of fields corresponding to precomputations.       
  dataSetInfo: Possible info about the dataset (useful when you store a model)
    globalOpt: A structure with the whole configuration file used for creating the model.
               You need to create the model with calling vargplvm_init.m and then use 
               the globalOpt to create the model if you want to actually use this
               global configuration -- this is not the creation way used in the tutorial
               but a simpler way is used. If you do want to use this global configuration code
               then check vargplvm_init.m.
     dynamics: The dynamics structure. This is the only structure that makes a Bayesian GP-LVM
               to turn into a VGPDS. If this structure is not present or is [], then the model
               and all vargplvm functions will treat the model as a BGPLVM, otherwise
               they will "understand" that it's actually a VGPDS. See below details.



%------ KERNEL (inside model.kern) ----------%

The kernel for the mapping. For the vargplvm model, we can use only kernels that make the
psi statistics tractable, ie:
Nonlinear kernels: 'rbfard2', 'rbf', or (recommended) 'rbfardjit' .
Linear kernels: 'linard2', 'lin'
For stability, we always use plus white plus (optionally) bias, by specifying a compound structure:
{'rbfard2','white','bias'};
rbfardjit is the exception, because it has the white and bias embedded.

You can even consider a linear and nonlinear kernel together:
{'rbfard2','linard2','white','bias'}.

In compound kernels, model.comp will have the compound structure and each cell is a kernel.

%- Example of an instantiation from the tutorial
>> model.kern
ans=
          type: 'rbfardjit' % It's like RBF + white noise + bias but more stable
inputDimension: 3
    transforms: [1x1 struct]
      variance: 0.7038
        jitter: 1.0000e-05
   inputScales: [0.0929 0.0865 2.3165e-04]
       nParams: 4
  isStationary: 1
        Kstore: []
         diagK: []



% Explanation of the fields:


          type: The type of kernel used for the mapping.
inputDimension: Dimensionality of inputs for k(x,x'), ie model.q usually.
    transforms: Transform the params of the kernel to e.g. force them to be positive
                (e.g. exponential transform)
      variance: The variance of the kernel
        jitter: This specific kernel (rbfardjit) has a constant jitter to make it numerically stable.
   inputScales: The ARD lengthscales associated with this kernel.
       nParams: Number of parameters
  isStationary: 
        Kstore: Don't worry about this field.
         diagK: Don't worry about this field.



%------ Variational distribution (inside model.vardist) ----------%

This implements q(X) = \prod_{n=1}^N GaussianDistribution(x_n | mu_n, S_n),
where mu_n is a vector and S_n is a diagonal matrix.

%- Example of an instantiation from the tutorial
>> model.vardist

ans = 

           type: 'vardist'
        vartype: 'gaussian'
        numData: 50
latentDimension: 3
        nParams: 300
     transforms: [1x1 struct]
          means: [50x3 double]
         covars: [50x3 double]


% Explanation of the fields:


           type: Leave as default (we don't have other types).
        vartype: Leave as default (we don't have other types).
        numData: Number of data, i.e. 
latentDimension: Dimensionality of q
        nParams: Redundant, keep because refactoring hasn't been done yet.
     transforms: You might want to constrain the parameters. usually we use this
                 to force the optimiser to find positive values for S_n
          means: A matrix where row j is the mean of q(x_j), ie mu_j
         covars: A matrix where row j is the variance of q(x_j),
                 ie the diagonal of S_j.



%------ DYNAMICS for VGPDS (inside model.dynamics) ----------%

This is the structure that turns a Bayesian GP-LVM into a model where the latent space
is constrained with a prior. If the prior is temporal (model.dynamics.constrainType = {'time'})
then we have VGPDS. But you might want to consider other constraints, e.g. you can create
a warped GP by arbitrary inputs Z to the prior for X, ie p(X | Z). For consistency we ALWAYS
refer to the mean of p(X) as model.dynamics.t, but it can obviously be whatever. For example,
we can even put model.dynamics.t = model.m to get an auto-encoder
(model.dynamics.constrainType = {'data'}) or we can create a Discriminative Bayesian GP-LVM
by constraining X with class labels (model.dynamics.constrainType = {'labels'}) or we can
even consider two or more constraints together
(model.dynamics.constrainType = {'labels','data'}), in which case the kernel for the latent space
(model.dynamics.kern) should be invcmpd (check that kernel for details).

If model.dynamics is not present as a field or is [], then the model is just a BGPLVM model.

%- Example of an instantiation from the tutorial

>> model.dynamics
ans = 

isMissingData: 0
            t: [50x1 double]
            X: [50x3 double]
            q: 3
            N: 50
         kern: [1x1 struct]
      vardist: [1x1 struct]
         type: 'vargpTimeDynamics'
 dynamicsType: 'regressive'
      nParams: 304
          seq: []
           Kt: [50x50 double]
   reoptimise: 1
learnVariance: 0
constrainType: {'time'}


% Explanation of the fields:

isMissingData: Don't worry about that field (leave default value).
            t: The input to the GP on the latent space. Can be time or something else. 
               It's a vector of size N or a matrix N x whatever.
            X: A copy of model.X.
            q: A copy of model.q
            N: A copy of model.n
         kern: The kernel for the top layer. This is DIFFERENT than model.kern. This kernel
               here is for the mapping model.t -> model.X. Since t is not a RV, there are 
               no Psi statistics involved, and this kernel can be whatever (see tutorial).
      vardist: This is NOT a copy of model.vardist. Due to the reparameterisation trick, the original
               mu and S of q(X) are stored in model.vardist (see above) and are reparametrised as
               bar(mu) and lambda, which are stored in model.dynamics.vardist. The model internally
               converts the reparemeterised parameters back to the original and vice versa whenever
               it's needed, automatically.
         type: We have only one type so far, 'vargpTimeDynamics'
 dynamicsType: We only have regressive dynamics so far.
      nParams: 
          seq: If [], then model.dynamics.t is treated as one sequence. Otherwise, if 
               e.g. seq = [s1 s2 s3] we will treat t(1:s1,:) as one sequence, t(s1+1:s1+s2, :) as
               another sequence etc, is each element in seq is the last index of model.dynamics.t
               for every sequence.
           Kt: The kernel matrix, computed by k_x(t,t)
   reoptimise: During test time, we have test latent points X*. Since everything is coupled,
               q(X*,X) is NOT factorised as q(X*) q(X) as in the BGPLVM case. therefore, by 
               making reoptimise field to be true, we reoptimise q(X) during test time (takes more
               time).
learnVariance: The dynamics kernel, model.dynamics.kern has the form sigma * k_x(..). Sometimes, we 
               might consider sigma to be redundant because of the lengthscale in k_x, and might want
               to fix it (eg. to 1).
constrainType: See discussion in the beginning of this section (Dynamics for VGPDS)











_____________________________ Demos and global configuration to run experiments

Check VGPDStutorial.m for a tutorial on VGPDS.
Check demosDynamics.m for more advanced demos on VGPDS.
Check demOilVargplvm4.m for a simple demonstration of Bayesian GP-LVM on
  the oil data.

Check MRDtutorial.m for a tutorial on MRD.

Check vargplvmEmbed2.m for the closest to using vargplvm and VGPDS as
   black boxes (perhaps checking VGPDStutorial.m will help) and
 and MRDEmbed.m and demSvargplvmGeneric.m. for the closest to using MRD
   as a black box (again, checking MRDtutorial.m is recommended).














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

A. Damianou, C. H. Ek, M. Titsias, N. Lawrence, Manifold Relevance
Determination, ICML 2012


A. Damianou, M. Titsias, N. Lawrence, Variational Gaussian process
dynamical systems, NIPS 2011


M. Titsias, N. Lawrence, Bayesian Gaussian process latent variable model, 
AISTATS 2010




___________________________________________________________________________
################ 5. Manifold Relevance Determination (MRD) ################
___________________________________________________________________________

Please consult svargplvm_README.txt


eof
