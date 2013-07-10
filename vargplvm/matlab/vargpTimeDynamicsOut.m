%%% this function is for test purposes, it's not used
function y = vargpTimeDynamicsOut(model, x);
    y = gpPosteriorMeanVar(model, x);
