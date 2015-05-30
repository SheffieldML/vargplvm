function scales = svargplvmShowScales(model, printPlots, modalities, dims, thresh, normaliseScales, customLeg, logScales)

% SVARGPLVMSHOWSCALES Show the scales of a svargplvmModel graphically
%
% SEEALSO : svargplvmFindSharedDims
%
% COPYRIGHT: Andreas C. Damianou, 2012, 2013
%
% VARGPLVM


if nargin < 2 || isempty(printPlots)
    printPlots = 2;
end

if nargin < 3 || isempty(modalities)
    modalities = 1:length(model.comp);
else
    printPlots = true;
end

if nargin < 4 || isempty(dims)
    dims = 1:model.q;
end

if nargin < 5 || isempty(thresh)
    thresh = -Inf;
end

if nargin < 6 || isempty(normaliseScales)
    normaliseScales = false;
end

if nargin < 7
    customLeg = {}; % Custom legend
end

if nargin < 8 || isempty(logScales)
    logScales = false;
end

totalModels = length(modalities);

if printPlots
    for i=1:length(model.comp)
        if strcmp(model.comp{i}.kern.type, 'rbfardjit')
            scales{i} = model.comp{i}.kern.inputScales(dims);
        else
            scales{i} = model.comp{i}.kern.comp{1}.inputScales(dims);
        end
    end
    if logScales
        origScales = scales;
        for k=1:totalModels
            scales{modalities(k)} = scales{modalities(k)}./max(scales{modalities(k)});
            scales{modalities(k)} = log(scales{modalities(k)});%logN(scales{modalities(k)},1.01);
            a=0; b=1; % Scale data between a and b
            X=scales{modalities(k)};
            X = (X-min(min(X))).*(b-a)./(max(max(X))-min(min(X))) + a;
            scales{modalities(k)}=X;
        end
    end
    if (totalModels == 2 || normaliseScales) && ~logScales
        origScales = scales;
        for k = 1:totalModels
            scales{modalities(k)}=scales{modalities(k)}./max(scales{modalities(k)});
        end
    end
    for k=1:totalModels
       scales{modalities(k)}(scales{modalities(k)}<thresh) = 0;
    end
    if totalModels == 2
        maxScales1 = max(scales{modalities(1)});
        maxScales2 = max(scales{modalities(2)});
        %scales{modalities(1)} = scales{modalities(1)}./maxScales1;
        %scales{modalities(2)} = scales{modalities(2)}./maxScales2;
        x=1:size(scales{modalities(1)},2); 
        K=0.5; 
        bar1=bar(x, scales{modalities(1)}, 'FaceColor', 'b', 'EdgeColor', 'b'); 
        set(bar1,'BarWidth',K); 
        hold on;
        bar2=bar(x, scales{modalities(2)}, 'FaceColor', 'r', 'EdgeColor', 'r');
        set(bar2,'BarWidth',K/2); 
        hold off; 
        if isempty(customLeg)
            legend(['scales1 (x ' num2str(maxScales1) ')'],['scales2 (x ' num2str(maxScales2) ')'])
        else
            legend(customLeg);
        end
        scales = origScales;
    else
        % Too many models. Plot one at a time
        %if totalModels > 10 && exist('hsvargplvmShowScales','file')
             %mm.layer{1} = model;
             %mm.H=1;
             %mm.layer{1}.M = totalModels;
             %hsvargplvmShowScales(mm);
        %else
            for i=1:totalModels
                if totalModels > 10
                    pause
                    %figure;
                elseif totalModels > 1
                    subplot(1,totalModels, i) ; 
                end
                bar(scales{modalities(i)}); title(['Modality ' num2str(modalities(i))])
            end
        %end
    end
else
    for i=1:length(model.comp)
        if strcmp(model.comp{i}.kern.type, 'rbfardjit')
            fprintf('# Scales of model %d: ', i)
            fprintf('%.4f  ',model.comp{i}.kern.inputScales(dims));
            fprintf('\n');
            scales{i} = model.comp{i}.kern.inputScales(dims);
        else
            fprintf('# Scales of model %d: ', i)
            fprintf('%.4f  ',model.comp{i}.kern.comp{1}.inputScales(dims));
            fprintf('\n');
            scales{i} = model.comp{i}.kern.comp{1}.inputScales(dims);
        end
        scales{i}(scales{i}<thresh) = 0;
    end
end


function y=logN(x,N)
y = log(x)/log(N);