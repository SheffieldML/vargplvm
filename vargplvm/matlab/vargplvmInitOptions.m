% Given a global options structure, create a local, equivalent one.

function options = vargplvmInitOptions(globalOpt, options)

if nargin < 2, options = []; end

if isempty(options)
    options = vargplvmOptions('dtcvar');
end

% See vargplvmOptions for what the following fields do
options.kern = globalOpt.mappingKern;
options.optimiser = globalOpt.optimiser;
options.initX = globalOpt.initX;
options.fixInducing = globalOpt.fixInd;
options.KLweight = globalOpt.KLweight;
options.initSNR = globalOpt.initSNR;
options.numActive = globalOpt.indPoints;