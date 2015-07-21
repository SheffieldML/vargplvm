function t_star = vgpdsPlotFit(model, Y, dims, limDf)

assert(isfield(model, 'dynamics') && ~isempty(model.dynamics));

if nargin < 2 || isempty(Y), Y = model.y; end
if nargin < 3 || isempty(dims), dims = 1:model.d; end
if nargin < 4 || isempty(limDf), limDf = 10; end

assert(length(dims) < 11); 

df = mean(diff(model.dynamics.t));
minT = min(model.dynamics.t) - limDf*df;
maxT = max(model.dynamics.t) + limDf*df;

t_star = linspace(minT, maxT, 300)';

% kk=3;t_star=[];
% for tt = 1:model.N-kk
% tmp=linspace(model.dynamics.t(tt),model.dynamics.t(tt+kk), 2*kk)';tmp=tmp(1:end);
% t_star = [t_star;  tmp];
% end
% t_star = unique(t_star);

[X, varX] = vargplvmPredictPoint(model.dynamics, t_star);
[Ypred, covYpred] = vargplvmPosteriorMeanVar(model, X, varX);

[p, ~] = numSubplots(length(dims) - 1);
numRows = p(1);
numColumns = p(2);

if numRows * numColumns < length(dims)
    [p, ~] = numSubplots(length(dims));
    numRows = p(1);
    numColumns = p(2);
end

if exist('tight_subplot', 'file') == 2
  ha = tight_subplot(numRows,numColumns,[.01 .03],[.1 .01],[.01 .01]);
end

for i=1:length(dims)
    d = dims(i);
    maxYdiff = abs(max(Y(:,d))) - abs(min(Y(:,d)));
    if exist('tight_subplot', 'file') == 2
        axes(ha(i));
    else
        subplot(numRows, numColumns,i);
    end
    plot(t_star, Ypred(:,d)); hold on
    plot(model.dynamics.t, Y(:,d), 'ro');
    legend('Fit', 'Training');
    covUp = Ypred(:,d)+sqrt(covYpred(:,d));
    covDown =  Ypred(:,d)-sqrt(covYpred(:,d));
    %covUp(covUp > Ypred(:,d) + maxYdiff) = NaN;
    %covDown(covDown < Ypred(:,d) - maxYdiff) = NaN;
    plot(t_star, covUp, ':');
    plot(t_star,covDown, ':');
    title(['d = ' num2str(d)]);
end


function [p,n]=numSubplots(n)
% Copyright (c) 2010, Rob Campbell
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.


% function [p,n]=numSubplots(n)
%
% Purpose
% Calculate how many rows and columns of sub-plots are needed to
% neatly display n subplots. 
%
% Inputs
% n - the desired number of subplots.     
%  
% Outputs
% p - a vector length 2 defining the number of rows and number of
%     columns required to show n plots.     
% [ n - the current number of subplots. This output is used only by
%       this function for a recursive call.]
%
%
%
% Example: neatly lay out 13 sub-plots
% >> p=numSubplots(13)
% p = 
%     3   5
% for i=1:13; subplot(p(1),p(2),i), pcolor(rand(10)), end 
%
%
% Rob Campbell - January 2010
   
    
while isprime(n) & n>4, 
    n=n+1;
end

p=factor(n);

if length(p)==1
    p=[1,p];
    return
end


while length(p)>2
    if length(p)>=4
        p(1)=p(1)*p(end-1);
        p(2)=p(2)*p(end);
        p(end-1:end)=[];
    else
        p(1)=p(1)*p(2);
        p(2)=[];
    end    
    p=sort(p);
end


%Reformat if the column/row ratio is too large: we want a roughly
%square design 
while p(2)/p(1)>2.5
    N=n+1;
    [p,n]=numSubplots(N); %Recursive!
end
