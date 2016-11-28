```
% SVARGPLVM_README
% Copyright: Andreas C. Damianou, 2012, 2013
% VARGPLVM
                                 _                 
                                | |                
  _____   ____ _ _ __ __ _ _ __ | |_   ___ __ ___  
 / __\ \ / / _` | '__/ _` | '_ \| \ \ / / '_ ` _ \ 
 \__ \\ V / (_| | | | (_| | |_) | |\ V /| | | | | |
 |___/ \_/ \__,_|_|  \__, | .__/|_| \_/ |_| |_| |_|
                      __/ | |                      
                     |___/|_|                      






                              R E A D M E

___________________________________________________________________________
```

## 1. GENERAL 
___________________________________________________________________________

This is the implementation of a Manifold Relevance Determination model, or
a "shared variational GP-LVM" (svargplvm). 
This model has been presented in the
Manifold Relevance Determination paper by A. Damianou, C. Ek, M. Titsias
and N. Lawrence, ICML 2012.

Similarly to the Bayesian GP-LVM, in the MRD the latent
space is marginalised out and the use of Automatic Relevance Determination
(ARD) covariance functions for the mapping from the latent to the observa-
tion space allows for automatic dimensionality detection. However, in
contrast to the Bayesian GP-LVM, here a different Gaussian Process (GP)
mapping is used per output modality, each with a different covariance
function and, hence, different set of ARD hyperparameters.

Notation: latent space X, output modalities Y^{(1)}, ..., Y^{(M)}

The MRD method consists of tools that allow for efficient initialisation
and training of such a model, so that after optimisation the different
sets of ARD hyperparamters define a "soft" segmentation for the latent
space, in the sense that each of the different modalities is allowed to
assign different importance to each dimension of the latent space, even
switching off completely some of those. If the ARD hyperparameters for
modality Y^{(1)} and Y^{(2)} are the same for dimension q, then we say that these
modalities share this dimension. Otherwise, dimension "q" can only be
active for one of the modalities ("private space") or for none (completely
switched off dimension). In general, even if dimension q is shared for,
e.g. modalties Y^{(1)} and Y^{(2)}, the corresponding weights 
w^{(1)}_1 and w^{(2)}_q can be dissimilar, hence the aforementioned "soft segmentation".

We might also choose one modality
per output dimension, i.e. we will have as outputs SEPARATE [y1 y2 y3...].
where subscript indexes dimensions.

! NOTE: This package contains a lot of files that are just serving the
purpose of utilities, ie most users will likely not need them. The core of
the necessary files is not as large. Check the demos and section 4 of this
README file.

! NOTE: This package heavily depends on the Bayesian GP-LVM (vargplvm soft-
ware package), you might want to spend some time getting familiar with
that one first.


___________________________________________________________________________
## 2. LATENT SPACE PRIORS 
___________________________________________________________________________

In the implementation, every output modality sees the same latent space.
It's the separate ARD weights that define the segmentation. For all sorts
of inference, these weights are automatically taken into account, so that
there is no need for any "postprocessing" step where we look at the weights
and assign latent spaces to the model.

The latent space can be constrained with priors (common, since as already
mentioned, the latent space is common). 

* The default prior is a standard normal, X ~ N(0, I), similar to the one
used in the Bayesian GP-LVM. 

* One can also use framework of VGPDS (Damianou et al. 2011) and constrain
the latent space with a temporal prior, i.e. the prior is now another GP
which warps its given input. This input can be 'time' (also given an asso-
ciated time vector), 'labels' (also given the training labels) or anything
else. The kind of constraints to use in the latent space is defined in a
field "dynamicsConstrainType" as a cell array (see svargplvm_init.m).


___________________________________________________________________________
## 3. WHAT CAN THIS MODEL DO? 
___________________________________________________________________________
Given a set of output spaces, Y^{(1)}, ..., Y_^{(M)}, the model comes up
with an (approximate) posterior q(X) ~ p(X | Y^{(1)}, ..., Y^{(M)}) and
associated ARD weights w^{(1)}, ..., w^{(1=M)}. Then:

a) the model can reveal the commonalities between the output spaces,
through the latent space. Just run the model and then visualise the ARD
weights (svargplvmShowScales(model)).

b) Exploring the latent space / correspondence problem:
the model can find the Nearest Neighbour similarities of one modality
relatively to another. For example, if X is segmented into
X = [X^{(1)} X^{(1,2)} X^{(2)}] (where X1 is private for modality 1, 
X2 for 2 and X12 is shared for both), we can give each training point y_n,
find the corresponding x_n and then from the shared dimensions of x_n (ie
the columns that correspond to the shared space) we can run the nearest
neighbour method to find those points of modality 2 that are the most
similar.

c) Similar to the above, but we give a TEST point in one modality and
produce a NOVEL output in the other modality. This is done by giving a test
point y^{(1)}_* in modality 1, for example, and finding the corresponding
q(x_*). Then, using this x_* (the mean of q(x_*)) we can generate a point
in the other modality, y^{(2)}_*. See section "Inference" in the paper to
see how this is done.

d) Classification.
Given a dataset (or multiple output modalities) and associated training
labels, we can treat the labels as an extra modality. Then, inference
can be carried on as in the case c) described above.

e) ...


* PREDICTIONS in detail
_______________________________________________________________________
A common prediction task would be to have the following scenario:
Assume we have 2 submodels in the svargplvm and after optimisation we end up
    with 2 private latent subspaces and 1 shared. E.g., the one model corresponds
    to a "figure" dataset, and the other model corresponds to the accompanying
    "pose" dataset. As regards predictions, we could give a point in the observed
    space of the one dataset (let's say, a pose y_i) and recover a point in the observed
    space of the other (a figure z_i). One way to do that is by first finding the
    corresponding latent points in the private space of the 1st model via the posterior
    q(X_star) ~ p(X_star|Y_star) -remember that q(X_star) has to be optimised variationally-.
    Now we have an X_* which is generated based on the private subspace of Y. To get an X which
    also takes into account the shared dimensions, we can take that X_* and find the NN of
    the latent points of the training data. Then, the prediction for Z will be given by the
    posterior.



___________________________________________________________________________
## 4. SOFTWARE DETAILS 
___________________________________________________________________________


To summarize, in the full model we can select:
   
* PARAMETERS for the actual MODEL:
    - The number of output modalities, M.
    - The prior distribution of the latent space, see section 2 above.
     - The mappings F^{(i)}, i=1:M between X and Y^{(i)} => 
                 kernel parameters and inducing points for each node
     - variational parameters for each node.


* OPTIMISATION:
    _______________________________________________________________________
    Given an initialised model (see below), the software runs a gradient-
    based optimiser which requires an objective function (see paper) and
    the associated derivatives to find a local optimum of the objective.
    The objective is a variational lower bound to the true evidence:
    log(Y^{(1)}, Y^{(2)}, ... , Y^{(M)}) and the parameters to optimise are:
        * \mu, \Sigma (the variational means and covariances for q(X)).
        * \tilde{X}^{(1)}, ..., \tilde{X}^{(M)} (the inducing points per modality)
        * \theta^{(1)}, ..., \theta^{(M)} (the kernel hyperparameters - including
            ARD weights- per modality)
        * \beta_1, ..., \beta_M (the inverse noise variance per modality, 
          coming from:
          y_m = f_m(X) + \epsilon_m, \epsilon_m ~ N(0, \beta_m^{-1}I)


* PARALLELISM:
    _______________________________________________________________________
    There is the option of parallelising the optimisation using the MATLAB
    parallel toolbox. The parallelism can be done in two ways (not recom-
    mended to do this for both at the same time):
    
    a) In the sub-model structures. If there are many modalities, then the
    bound / gradient computations for each can be done separately.
    To activate this, pass a field model.comp{i}.parallel = 1.
    
    b) In the variational distribution. If there are many datapoints, then
    the computation of the Psi statistics (see Bayesian GP-LVM and the
    VARGPLVM toolbox) can be done parallely w.r.t datapoints. This needs
    model.vardist.parallel = 1.

    Check svargplvmPropagateField.m to see how flags can be set/unset for
    all submodels at once.

    So, paralellism can be achieved in two levels: with respect to models, and with
    respect to N in each model (this is done in the rbfard2VardistPsi2Compute
    and rbfard2VardistPsi2Gradients). 

    We should perform some benchmarks and see which are the optimal settings,
    but one would expect that it is not optimal to use both ways to parallelize
    the program. 
    A reasonable idea is to use the model-wise parallelism when we have a large
    number of models (e.g. close to the number of cores/workers to be used),
    otherwise (and if N is relatively big) we can use the point-wise paral


* INITIALISATION:
    _______________________________________________________________________
    The model is prone to "bad" local optima. Therefore, initialisation is
    important. Most of the initialisation options can be explored by
    looking at the initialisation function "svargplvm_init.m" and studying
    the demos / tutorial. 
    
    The initialisation is done as follows:
    Every demo is a script, and requires a list of variables that define
    the initial conditions of the demo. This list is created by the
    script "svargplvm_init.m". This script returns a struct "globalOpt"
    which holds all initialisation options. To do that, svargplvm_init.m:
        checks for every field 'F', if there is a variable 'F' in the
        workspace. If yes, it copies globalOpt.F the value found. Otherwise
        it copies a default value.

    In brief, the things to take care of during initialisation, are:
        * The way in which the latent space is initialised. See
        the field "initial_X" in svargplvm_init.m
        * The initial \mu, \Sigma, inducing points, ARD weights. These
        initialisations should be taken care of automatically in the
        software.
        * If initialisation is bad, the model may find the trivial local
        optimum where everything is explained by noise (ie \beta == var(Y)).
        To avoid this, you can optimise the model for a few iterations with
        the \beta parameter being fixed. Then, this parameter can be
        "released" to the optimiser. This is taken care of automatically in
        the software, by specifying the number of iterations for
        initialising the variational distribution in that way (field
        "initVardistIters"). Then, field "itNo" specifies the number of
        "normal" optimiser iterations to run. See next bullet for SNR.
        

    * Signal to Noise Ratio (SNR): The last bullet above, practically refers
      to fixing the Signal to Noise Ratio of the model for a few iterations,
      forcing it to learn less noise in the first few iterations.
      The initial SNR for the constraint iterations (where a fixed percentage
      of noise is learned by keeping beta fixed) is 100, but you can set it
      yourself in the same way that all initial fields are set, (see
      svargplvm_init.m and demos).


* MODEL STRUCTURE AND FIELDS:
    _______________________________________________________________________
    The "model" structure holds all information for the model.
    There is a substructure model.comp which is a cell array, one cell per
    modality. 

    * model.comp: Each cell, is a sub-model of type vargplvm (see the VARGPLVM
    package which implements the Bayesian GP-LVM). However, these sub-models
    share the same latent space. In brief, this submodel contains the basic
    elements of a vargplvm, such as "m", (the centered version of the data),
    beta, X_u (the inducing points), the structure "kern" which implements
    the kernel, K_uu which is the covariance matrix evaluated on the ind.
    points, and many other precomputations used for computing the objective
    function (see below) and the dericatives.

    * model.vardist: This latent space is stored in
    model.vardist, where model.vardist.means is \mu (rows are points, columns
    dimensions 1,..., Q) and model.vardist.covars is \Sigma. Notice that
    \Sigma is stores the diagonals of all individual \Sigma_q. We don't
    need to store the full matrix \Sigma_q, as this is either diagonal
    anyway (if the latent space prior is standard normal), or the off-diagonal
    elements are only needed for the warping GP and are used only in the
    dynamics sub-field which implements this GP (see below).

    * model.dynamics: If a non-standard normal prior is used for X, then
    this field will be nonempty and will implement the dynamics (or any
    sort of prior) as described in the VGPDS paper (Damianou et al. 2011).
    This structure will also contain model.dynamics.Kt, ie the covariance
    matrix built on the prior component which takes known inputs which
    are denoted by "t" but can be something other than timestampts as well.

    * Details:
    The model is basically a collection of variational GPLVM models.
    e.g.
```
    model = 
    
            comp: {[1x1 struct]  [1x1 struct]} % comp{i} is the vargplvm model i.
               N: 100	% number of datapoints in each dataset (has to be the same)
               q: 6     % number of latent dimensions
         vardist: [1x1 struct]   % shared var. distr q(X)
               X: [100x6 double] % shared latent points (means of q(X))
       numModels: 2		% number of vargplvm models
            type: 'svargplvm'
    dataSetNames: {'fols_cos'  'fols_sin'}
      initLatent: 'pca'
    experimentNo: 404
        dataType: 'fols'
     initVardist: 1
```
	 
Each of the vargplvm submodels keeps an updated version of the shared elements (e.g.
q(X)), because there is abstraction, i.e. the vargplvm models don't "know" that they
are optimised within a svargplvm rather than within a vargplvm framework.

!! Check vargplvm_README.txt to see the exact structure of the "model"
    used there. That will be helpful, since MRD is essentially a model
    where every model.comp{i} element is a modality represented by a whole
    vargplvm (Bayesian GP-LVM) model.


* FUNCTIONS, FILES, DEPENDENCIES:
    _______________________________________________________________________
    * Objective function: svargplvmLogLikelihood.m
    * Associated derivatives: svargplvmLogLikeGradients.m
    * The optimiser function sends and receives parameters as a vector,
      but the above functions operate by directly getting the parameters
      from the model structure. Therefore, wrapper functions
      svargplvmExtractParam.m and svargplvmExpandParam.m are used for
      extracting and expanding the parameter vector from and into the model.

    ...

* UTILITIES:
    _______________________________________________________________________
    There are a few files that are not vital for the method itself, but are
    useful utilities:
     - svargplvmPropagateField.m
     - svargplvmPruneModel.m / svargplvmRestorePrunedModel.m
     - svargplvmShowScales.m
     - svargplvmSNR.m
     - svargplvmEqualizeScales.m
     - svargplvmFindSharedDims.m
     - ...

___________________________________________________________________________
## 4. DEMOS / TUTORIALS 
___________________________________________________________________________
    
* 4.1 Demos / Tutorials:
    _______________________________________________________________________
    There are plenty of demos, in the folder demos/. The demos start with
    the prefix "dem". Some of the demos are more well organised, some are
    mostly there for debugging purporses / depricated and might have bugs
    and inconsistencies. The best way to start with the software is to
    run the tutorial MRDtutorial.m, which calls the various demos. Then
    you can experiment with the demos by tweaking them, etc.

* 4.2 Can I use MRD as a "black box"?
    _______________________________________________________________________
    Possibly, by taking one of the demos and replacing the data structure
    with your data structure and fixining the required numbers (e.g. if
    your data is 5-dimensional, it makes sense to initialise X with less than
    5 dimensions). 
    However, it's very likely that attempting to use MRD as a black box may
    result in unsatisfactory or even bad results, since:
        a) some basic understanding is needed to avoid inconsistencies
           like the one mentioned above, e.g. if D is 10, X cannot be > 10.
           Or, for example, when attempting to model a highly non-linear
           dataset with a linear mappping kernel.
        b) to initialise the model in a reasonable way. This doesn't mean
           that you are required to do all sorts of hacks, but for challenging
           datasets it should help.
        c) there are many options to initialise the demos. By default, the most
           common are used, but there is the possibility that one option
           is better than the other for specific datasets. For example, the
           choice of mapping / warping kernel matters.
    



```
_____________________________ ChangeLog
v. 1.0: First release. Last update: 1/4/2013



_____________________________ TODO
The variational distribution q(X) should be changed so that it's included
only once in the model structure, instead of existing in every sub-model
in model.comp.
```

