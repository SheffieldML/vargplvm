function [gVarmeans gVarcovs gDynKern] = ...
    vargpTimeDynamicsVarPriorGradients(model, gPsi1, gPsi2, gPsi0)

% VARGPTIMEDYNAMICSVARPRIORGRADIENTS This function does all the work for calculating the
% derivatives for the variational parameters (the new ones, lambda, mu_bar)
% and for the theta_t parameters. The order of the calculations is such
% that everything is already computed when it's needed without
% recalculating things (e.g. Sq's are only calculated once). 
%
% THIS FUNCTION IS NOT USED AND WILL BE REMOVED FROM FUTURE RELEASES!!!!!
% Notes: model.dynamics.vardist is used to hold lambda_q's and mu_bar's.
% The original mu's and Sq's are not stored anywhere, just computed
% locally, (Sq and mu_q). In that sense, model.vardist is probably not
% necessary. 
%
% Most of the things are doable, only where I've put "TODO!!!" there is a
% little problem that must be taken care of.
%
% VARGPLVM

Kt = kernCompute(model.kern, model.t);
Lt = jitChol(Kt)';
invLt = Lt \ eye(model.N);
invLtT = invLt';
invKt = invLtT * invLt;

sumTr_q = 0;
sumgKt_q = 0;


for q=1:model.q
     % Calculate S_q
     Sqinv= invKt + diag(model.vardist.covars(q));
     Sq = Sqinv \ eye(model.N); %better with cholesky
     
     %------ Calculate gKL (report eq. 132)
     % Precomputations (A0_q and a1_q)
     Ls_q = jitChol(Sq)';
     model.mu{q} = Kt * model.vardist.means(q);     
     a1_q = invLt * model.mu{q};
     A0_q = invLt * Ls_q;
     
     % Calculate the first term of the trace
     tempTerm = invLtT * A0_q;
     tempTerm = - sum(sum(tempTerm .* tempTerm));
     sumgKt_q = sumgKt_q + tempTerm;
     
     % ... second term
     tempTerm = invLtT * a1_q;
     tempTerm = sum(sum(tempTerm .* tempTerm));
     sumgKt_q = sumgKt_q + tempTerm;
     % ... third term
     sumgKt_q = sumgKt_q + invKt;
     
     % kernGradient will now trace sumgKt_q after multiplying with the
     % derivative of Kt w.r.t the parameter theta_t. Add this to the big
     % summation over q.
     sumTr_q = sumTr_q + kernGradient(model.kern, model.t, sumgKt_q);
     sumTr_q = -0.5 * sumTr_q;
     
     %------ Calculate the deriv. of the likelihood term w.r.t diag(S_q).
     % the following needs changing (see "TODO")
     [ gVarmeans1, gVarcovs1] = kernVardistPsi1Gradient(model.dynamics.vardist, gPsi1');
     [ gVarmeans2, gVarcovs2] = kernVardistPsi2Gradient(model.dynamics.vardist, gPsi2);
     [ gVarmeans0, gVarcovs0] = kernVardistPsi0Gradient(model.dynamics.vardist, gPsi0);
      % TODO!!!!!! Need modified functions which operate only on sq, mu_q!!! 
     
     %------ Calculate the deriv. of the whole bound w.r.t lambda_q 
     %------(report eq. (98)
     % (S1 .* Sq) * (gVarcovsLik_q - * lambda_q).
     % where lambda_q = model.dynamics.vardist.covars(:,q)
     
     
     %------ Calculate the deriv. of the whole bound w.r.t mu_bar_q
     %------- (report eq. 90)
     % Kt * (gVarmeansLik_q - mu_bar_q)
     % where mu_bar_q = model.dynamics.vardist.means(:,q)
     
     %------ Calculate the deriv. of the whole bound w.r.t a kernel param. of
     %------ Kt (amending the derivative of KL).  (report eq.103)
        % TODO!!!
        % I miss a sum from the formulas??? It shouldn't be just ONE
        % mu_q...
        % Otherwise, I can use gKernKL, gVarmeansLik_q,  gVarcovsLik_q, Sq,
        % Kt, all already computed. Finally, I'll also need the deriv of Kt
        % w.r.t theta_q (it is already implemented).
end
gKernKL = sumTr_q;
