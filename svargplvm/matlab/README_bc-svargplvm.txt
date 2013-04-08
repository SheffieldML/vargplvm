% This runs the demo for the svargplvm model which is as follows:
% theta_x -> X -> {Y, Z}, where Y,Z are the two output modalities
% (one can put Y as the data and Z as the corresponding labels in
% 1-of-K encoding with 1s and 0s or with 1s and -1s) and theta_x
% is whatever constrains X. It can be either nothing, ie standard
% normal p(X) (no "dynamics"), or 'data' (back-constrained model)
% or 'labels' (if provided) or 'time' (if provided).
% Check the bc-vargplvm package for how to add multiple constraints
% using "inverse of sum of inverses" kernels.

% Note: this demo needs the bc-vargplvm package added. It is similar
% to the bc-vargplvm package, but it allows multiple output modalities.

% This is a general demo, change the options here according to your
% needs / specific dataset you want to use.

   dataSetName = 'mulan_emotions';
% For the following options, you can also check the comments in
% svargplvm_init which explain all the possible options

   % The mapping kernel
   baseKern = {{'linard2','white'}, {'linard2','white'}}; 
   indPoints = 85;
   initX = 'pca';
   % Initialise by concatenating the output modalities and then
   % doing PCA (check svargplvm_init about how to do it differently)
   initial_X = 'concatenated';
   latentDim = 10;
   % Initial SNR
   initSNR = 100;
   % Optimisation iterations with beta being fixed
   initVardistIters = 500;
   % Sometimes doing some iterations (eg 500) and then restarting the optimiser
   % and doing some more (eg 200) is better than doing all of them (eg 700)
   % altogether, as it helps avoiding local minima. If left blank ([]) we don't
   % do this optimisation at all, ie we only do optimisation for "initVardistIters"
   % with SNR fixed
   itNo = [500 200];
   % Scale data to have variance one
   scale2var1 = 1;  
   % The dynamical model is theta_x -> X -> Y.
   % The following field says what theta_x is. It can be
   % 'time', 'labels', 'data' or {} (no "dynamics", just a standard normal
   % on X).
   dynamicsConstrainType = {'labels'}; % dynIsed = 1
   % The dynamical kernel to use
   dynamicKern = {'lin', 'white', 'bias'};
   % This has to be set to a value so that model.vardist.covars is around 0.3
   vardistCovarsMult = 1;
   doPredictions = true;   
 
   demClassificationGeneral


