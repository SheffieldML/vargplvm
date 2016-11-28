% DEM_MRD_SIMPLE
% Simpler way of running MRD. This demo is embedding the modalities given
% in Yall into a latent space X which is returned. The MRD model is also
% returned. Parameters for the MRD are given as arguments in the call for
% MRDEmbed (See function for arguments).
%
% See MRDtutorial.m for a more complete version of this demo.
%
% COPYRIGHT: Andreas Damianou 2016

% Get some toy data
Yall = util_createMRD_toy_data();
%kernel to use
baseKern = {{'linard2','white', 'bias'},{'linard2','white', 'bias'}}; 
% Learn MRD with some arguments.
[X,model] = MRDEmbed(Yall, 'initial_X','separately','latentDimPerModel',4,'baseKern',baseKern,'initVardistIters',200,'itNo',200);
