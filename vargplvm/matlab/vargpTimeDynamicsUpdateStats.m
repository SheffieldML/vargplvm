function model = vargpTimeDynamicsUpdateStats(model)
% VARGPTIMEDYNAMICSUPDATESTATS Supplementary update stats for when the model contains
% dynamics
% FORMAT
% DESC When the model contains dynamics, this function first updates the
% dynamics structure with the relevant precomputations (including the
% calculation of the temporal prior's kernel Kt) and then it
% calculates the model's original variational means and covariances. This
% is essential because due to the reparametrization strategy followed (see
% report) the optimiser and the likelihood-gradient functions operate on
% the newly introduced parameters which are, however, a function of the
% original ones. The original ones must, thus, be computed in each cycle.
% Note that the VARDIST field of the dynamics structure contains the newly
% introduced parameters and the VARDIST field of the model structure
% contains the original parameters.
%
% The code supports the case when the model learns from multiple sequences. 
% In that case, instead of calculating Kt as usual (given the vector t) it
% forces it to be block diagonal where each block corresponds to a specific
% sequence. This assumes, though, that correlations between sequences are
% not modelled.
% 
% COPYRIGHT : Andreas C. Damianou, 2010-2011
% COPYRIGHT : Michalis K. Titsias, 2010-2011
% COPYRIGHT : Neil D. Lawrence, 2010-2011

% 


% VARGPLVM




model.dynamics.X = model.vardist.means;

%if isfield(model.dynamics, 'fixedKt') & ~isempty(model.dynamics.fixedKt)
if isfield(model.dynamics, 'fixedKt')  & (model.dynamics.fixedKt == 1)
    % do nothing assume that Kt is correct    
%   model.dynamics.Kt = model.dynamics.fixedKt;
else
    if isfield(model.dynamics, 'seq') & ~isempty(model.dynamics.seq)
        % Create a block diagonal Kt matrix according to indices. If seq == N
        % then this is equivalent to creating Kt out of the whole t vector.
        model.dynamics.Kt = zeros(model.dynamics.N);
        tindexPrev = 1;
        for i=1:length(model.dynamics.seq)
            tindexCur = model.dynamics.seq(i);
            tseq = model.dynamics.t(tindexPrev:tindexCur);
            Kt_i = kernCompute(model.dynamics.kern, tseq);
            model.dynamics.Kt(tindexPrev:tindexCur, tindexPrev:tindexCur) = Kt_i;
            tindexPrev = model.dynamics.seq(i)+1;
        end
    else
        model.dynamics.Kt = kernCompute(model.dynamics.kern, model.dynamics.t);
    end
end


%%%% Update the base model with the original variational parameters.
if isfield(model.dynamics, 'seq') & ~isempty(model.dynamics.seq)
    [model.vardist.means model.vardist.covars] = vargpTimeDynamicsPriorKernGradSeq(model.dynamics);
else
    [model.vardist.means model.vardist.covars] = vargpTimeDynamicsPriorKernGrad(model.dynamics);
end



function [muqOrig SqOrig] = vargpTimeDynamicsPriorKernGradSeq(dynModel)
%fprintf(1,'vargpTimeDynamicsVarPriorKernGrad2..\n');%%% DEBUG


muqOrig=zeros(dynModel.N, dynModel.q);
SqOrig = zeros(dynModel.N, dynModel.q); % memory preallocation

seq = dynModel.seq;
for q=1:dynModel.q
    LambdaH_q = dynModel.vardist.covars(:,q).^0.5;
    
    seqStart=1;
    for i=1:length(dynModel.seq)
        seqEnd = seq(i);
             %   size(dynModel.Kt) %%%%

        %dynModel.Kt(seqStart:seqEnd,seqStart:seqEnd); %%%!!
        Bt_q = eye(seqEnd-seqStart+1) + LambdaH_q(seqStart:seqEnd,1)*LambdaH_q(seqStart:seqEnd,1)'.*dynModel.Kt(seqStart:seqEnd,seqStart:seqEnd);
        Lbt_q = jitChol(Bt_q)';
        G = repmat(LambdaH_q(seqStart:seqEnd,1), 1, seqEnd-seqStart+1).*dynModel.Kt(seqStart:seqEnd,seqStart:seqEnd);
        G = Lbt_q \ G;
        SqOrig(seqStart:seqEnd,q) = diag(dynModel.Kt(seqStart:seqEnd,seqStart:seqEnd)) - sum(G.*G,1)';
        muqOrig(seqStart:seqEnd,q) = dynModel.Kt(seqStart:seqEnd,seqStart:seqEnd) * dynModel.vardist.means(seqStart:seqEnd,q);

        seqStart = seqEnd+1;
    end
end



% Original
function [muqOrig SqOrig] = vargpTimeDynamicsPriorKernGrad(dynModel)
%fprintf(1,'vargpTimeDynamicsVarPriorKernGrad2..\n');%%% DEBUG

muqOrig = dynModel.Kt * dynModel.vardist.means;

SqOrig = zeros(dynModel.N, dynModel.q); % memory preallocation

for q=1:dynModel.q
    LambdaH_q = dynModel.vardist.covars(:,q).^0.5;
    Bt_q = eye(dynModel.N) + LambdaH_q*LambdaH_q'.*dynModel.Kt;
    
    Lbt_q = jitChol(Bt_q)';

    G = repmat(LambdaH_q, 1, dynModel.N).*dynModel.Kt;
    G = Lbt_q \ G; 
    % Numerical stable way for diag(inv(inv(dynModel.Kt +
    % diag(dynModel.vardist.covars(:,q))));
    SqOrig(:,q) = diag(dynModel.Kt) - sum(G.*G,1)';
    
end


