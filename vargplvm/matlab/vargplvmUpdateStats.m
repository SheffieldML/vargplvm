function model = vargplvmUpdateStats(model, X_u)
% VARGPLVMUPDATESTATS Update stats for VARGPLVM model.
%
% COPYRIGHT : Michalis K. Titsias, 2009-2011
% COPYRIGHT : Neil D. Lawrence, 2009-2011
% COPYRIGHT : Andreas C. Damianou, 2010-2011
% 
% SEEALSO : vargplvmOptimise, vargplvmExpandParam

% VARGPLVM
  
jitter = 1e-6;
%model.jitter = 1e-6;


% %%% Precomputations for the KL term %%%
%%% Maybe we should add something like model.dynamics.X =
%%% model.dynamics.vardist.means (if visualisation requires that).
 if isfield(model, 'dynamics') && ~isempty(model.dynamics)
     % This actually only operates on model.dynamics. The only reason we
     % send the whole model is that we need to get model.vardist.means.
     model = vargplvmDynamicsUpdateStats(model);
 end



%%% Precomputations for (the likelihood term of) the bound %%%

model.K_uu = kernCompute(model.kern, X_u);

% Always add jitter (so that the inducing variables are "jitter" function variables)
% and the above value represents the minimum jitter value
% Putting jitter always ("if" in comments) is like having a second
% whiteVariance in the kernel which is constant.
%if (~isfield(model.kern, 'whiteVariance')) | model.kern.whiteVariance < jitter
   % There is no white noise term so add some jitter.
if ~strcmp(model.kern.type, 'rbfardjit')
       model.K_uu = model.K_uu ...
                         + sparseDiag(repmat(jitter, size(model.K_uu, 1), 1));
end
%end

model.Psi0 = kernVardistPsi0Compute(model.kern, model.vardist);
model.Psi1 = kernVardistPsi1Compute(model.kern, model.vardist, X_u);
[model.Psi2, AS] = kernVardistPsi2Compute(model.kern, model.vardist, X_u);

%K_uu_jit = model.K_uu + model.jitter*eye(model.k);
%model.Lm = chol(K_uu_jit, 'lower');          
  
                                      % M is model.k
%model.Lm = chol(model.K_uu, 'lower'); 
model.Lm = jitChol(model.K_uu)';      % M x M: L_m (lower triangular)   ---- O(m^3)
model.invLm = model.Lm\eye(model.k);  % M x M: L_m^{-1}                 ---- O(m^3)
model.invLmT = model.invLm'; % L_m^{-T}
model.C = model.invLm * model.Psi2 * model.invLmT;
model.TrC = sum(diag(model.C)); % Tr(C)
% Matrix At replaces the matrix A of the old implementation; At is more stable
% since it has a much smaller condition number than A=sigma^2 K_uu + Psi2
model.At = (1/model.beta) * eye(size(model.C,1)) + model.C; % At = beta^{-1} I + C

model.Lat = jitChol(model.At)'; % UNCOMMENT THIS IF YOU ARE NOT DOING THE DEBUG BELOW

%try %%%% DEBUG1
% model.Lat = jitChol(model.At)';
%catch e %%% DEBUG
%    model_At = model.At;
%    save('modellAt', 'model_At');
%    e.throw %%% DEBUG
%end %%%% DEBUG


% DEBUG2
% warning('') %%% DEBUG
% model.Lat = jitChol(model.At)';
% tmpWarn = lastwarn; %%% DEBUG
% if length(tmpWarn) > 41 && strcmp(tmpWarn(1:42), 'Matrix is not positive definite in jitChol') %%% DEBUG
%     % do stuff...
%     %%%model.id
%     warning('') % This is a good place to put a stop for the debugger...
% end %%% DEBUG



model.invLat = model.Lat\eye(size(model.Lat,1));  
model.invLatT = model.invLat';
model.logDetAt = 2*(sum(log(diag(model.Lat)))); % log |At|

model.P1 = model.invLat * model.invLm; % M x M

% First multiply the two last factors; so, the large N is only involved
% once in the calculations (P1: MxM, Psi1':MxN, Y: NxD)
model.P = model.P1 * (model.Psi1' * model.m);

% Needed for both, the bound's and the derivs. calculations.
model.TrPP = sum(sum(model.P .* model.P));


%%% Precomputations for the derivatives (of the likelihood term) of the bound %%%

%model.B = model.invLmT * model.invLatT * model.P; %next line is better
model.B = model.P1' * model.P;
model.invK_uu = model.invLmT * model.invLm;
Tb = (1/model.beta) * model.d * (model.P1' * model.P1);
	Tb = Tb + (model.B * model.B');
model.T1 = model.d * model.invK_uu - Tb;


model.X = model.vardist.means;

%{
%%%%-------------------------- TEMP (for analysing how ind. pts are optimised) !!!!!!!!!
d=dist2(model.X_u, model.vardist.means);
[C,I]=min(d,[],2);
dNorm=dist2(model.X_u, model.vardist.means)/norm(model.X_u)^2;
[Cnorm,Inorm]=min(dNorm,[],2);

if Inorm-I ~= 0
    fprintf(1, 'Something is wrong with the distances! (Update Stats)');
end

try 
    load TEMPInducingMeansDist
catch exception
    TEMPInducingMeansDist=[];
end
try 
    load TEMPInducingIndices
catch exception
    TEMPInducingIndices=[];
end
try 
    load TEMPInducingMeansDistNorm
catch exception
    TEMPInducingMeansDistNorm=[];
end
TEMPInducingMeansDist = [TEMPInducingMeansDist sum(C)];
save 'TEMPInducingMeansDist.mat' 'TEMPInducingMeansDist';
%fprintf(1,'# Total X_u - vardistMeans = %d\n',sum(C));
TEMPInducingIndices = [TEMPInducingIndices I];
save 'TEMPInducingIndices.mat' 'TEMPInducingIndices';

TEMPInducingMeansDistNorm = [TEMPInducingMeansDistNorm sum(Cnorm)];
save 'TEMPInducingMeansDistNorm.mat' 'TEMPInducingMeansDistNorm';
%}

