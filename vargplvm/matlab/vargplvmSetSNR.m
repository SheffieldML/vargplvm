function model = vargplvmSetSNR(model, SNR, displ)
if nargin < 3 || isempty(displ), displ = true; end

assert(nargin > 1)

if isfield(model, 'mOrig') && ~isempty(model.mOrig)
    varY = var(model.mOrig(:));
else
    varY = var(model.m(:));
end

beta = model.beta;
model.beta = SNR/varY;

if displ
    fprintf('     %f  (varY=%f, 1/beta=%f)\n',  SNR, varY, 1/beta)
end