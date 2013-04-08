function x_cur = svargplvmPredictX(x_star, X, kNN, sharedDims, options)

% SHEFFIELDMLPREDICTX Given a test latent point found by a corresponding output point in one modality, this function
% finds the latent point to be used for generating outputs in the other modality.
%
% SEEALSO : demClassification, demClassification3, demClassificationGeneral, svargplvmPredictions2, svargplvmPredictions3
%
% COPYRIGHT : Andreas C. Damianou, 2011

% SHEFFIELDML

switch options
    case 1
        x_cur = x_star;
    case 2
        x_cur = X(kNN,:);
    case 3
        x_cur = X(kNN,:);
        x_cur(sharedDims) = x_star(sharedDims);
        fprintf('--'); %%
end

