function pb = vargpTimeDynamicsVarPriorBound(dynModel)
% VARGPTIMEDYNAMICSVARPRIORBOUND Computes the term of the variational
% lower bound which corresponds to the dynamics part, where the prior is
% incolved.
% FORMAT
% DESC It takes a dynamic model structure and calculates the term of the
% bound that only corresponds to the dynamics part, i.e. to the KL term
% which involves a temporal prior. The bound is a lower
% Jensens bound for variational approximation. The code supports
% calculations for when there are multiple sequences to learn from.
%
% ARG dynModel : the dynamic model structure for
% which the corresponding part of the bound is to be computed.
% RETURN pb : the term of the variational bound which corresponds to the 
% dynamics structure. This term comes from a KL term.
% 
% SEEALSO : vargplvmLogLikelihood, modelVarPriorBound
%
% COPYRIGHT : Michalis K. Titsias, 2010-2011
% COPYRIGHT : Andreas C. Damianou, 2010-2011
% COPYRIGHT : Neil D. Lawrence, 2010-2011

% VARGPLVM


if isfield(dynModel, 'seq') & ~isempty(dynModel.seq)
    pb = vargpTimeDynamicsVarPriorBoundSeq(dynModel);
else
    pb = vargpTimeDynamicsVarPriorBound1(dynModel);
end


% The following function is only called if the model learns from multiple
% independent sequences. In that case, the code can be optimised because Kt
% (as well as the matrices built using that, such as Bt)
% is block diagonal and we can perform individual operations in each block
% rather than on the whole Kt matrix directly. 
function pb = vargpTimeDynamicsVarPriorBoundSeq(dynModel)

sumMu = zeros(dynModel.N, dynModel.N);

% Initialize the sums for each block (sequence)
seq = dynModel.seq;
seqStart=1;
for i=1:length(dynModel.seq)
    seqEnd = seq(i);
    sumBhat{i} = zeros(seqEnd-seqStart+1, seqEnd-seqStart+1);
    sumLogDetB{i}=0;
    seqStart = seqEnd+1;
end

% The following code implements the equations presented in the report,
% which are derived using the Matrix Inversion Lemma. In that way, the Kt
% matrix is never inverted, but instead, the more stable Bt matrix is.
for q=1:dynModel.q

    LambdaH_q = dynModel.vardist.covars(:,q).^0.5;

    seqStart=1;
    for i=1:length(dynModel.seq)
        seqEnd = seq(i);
        
        Bt_q = eye(seqEnd-seqStart+1) + LambdaH_q(seqStart:seqEnd,1)*LambdaH_q(seqStart:seqEnd,1)'.*dynModel.Kt(seqStart:seqEnd,seqStart:seqEnd);
        Lbt_q = jitChol(Bt_q)';
        Bhat_q = Lbt_q \ diag(LambdaH_q(seqStart:seqEnd,1));
        Bhat_q = Bhat_q'*Bhat_q;

        %sumBhat = sumBhat + Bhat_q;
        sumBhat{i} = sumBhat{i} + Bhat_q;
        sumLogDetB{i} = sumLogDetB{i} + 2*sum(log(diag(Lbt_q)));

        seqStart = seqEnd+1;
    end
    sumMu = sumMu + (dynModel.vardist.means(:,q) * dynModel.vardist.means(:,q)');
end

sumLogDetBAll=0;
tmp=eye(dynModel.N);
seqStart=1;
for i=1:length(seq)
    sumLogDetBAll=sumLogDetBAll + sumLogDetB{i};
      seqEnd = seq(i);
    tmp(seqStart:seqEnd, seqStart:seqEnd)=sumBhat{i};
    seqStart = seqEnd+1;
end

sumBhat = tmp;
sumLogDetB = sumLogDetBAll;
pb = 0.5*(sum(sum(sumBhat.*dynModel.Kt)) - sum(sum(dynModel.Kt .* sumMu)) - sumLogDetB ); 


% This function is equivalent to the above when no sequencing is used. All
% of the operations are performed in full matrices.
function pb = vargpTimeDynamicsVarPriorBound1(dynModel)


sumBhat = zeros(dynModel.N, dynModel.N);
sumMu = zeros(dynModel.N, dynModel.N);
sumLogDetB = 0;

for q=1:dynModel.q

    LambdaH_q = dynModel.vardist.covars(:,q).^0.5;

    Bt_q = eye(dynModel.N) + (LambdaH_q*LambdaH_q').*dynModel.Kt;
   
    Lbt_q = jitChol(Bt_q)'; 
    

    Bhat_q = Lbt_q \ diag(LambdaH_q);
    Bhat_q = Bhat_q'*Bhat_q; 
    sumBhat = sumBhat + Bhat_q; 

    sumMu = sumMu + (dynModel.vardist.means(:,q) * dynModel.vardist.means(:,q)');   
    sumLogDetB = sumLogDetB + 2*sum(log(diag(Lbt_q)));
end

% Return the value of the bound.
pb = 0.5*(sum(sum(sumBhat.*dynModel.Kt)) - sum(sum(dynModel.Kt .* sumMu)) - sumLogDetB ); 
      
