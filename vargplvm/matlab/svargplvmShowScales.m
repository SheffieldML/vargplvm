function scales = svargplvmShowScales(model, printPlots)

% SVARGPLVMSHOWSCALES Show the scales of a svargplvmModel graphically
%
% SEEALSO : svargplvmFindSharedDims
%
% COPYRIGHT: Andreas C. Damianou, 2012, 2013
%
% VARGPLVM

totalModels = length(model.comp);

if nargin < 2
    printPlots = true;
end

if printPlots
    for i=1:totalModels
        if strcmp(model.comp{i}.kern.type, 'rbfardjit')
            scales{i} = model.comp{i}.kern.inputScales;
        else
            scales{i} = model.comp{i}.kern.comp{1}.inputScales;
        end
    end
    if totalModels == 2
        origScales = scales;
        maxScales1 = max(scales{1});
        maxScales2 = max(scales{2});
        scales{1} = scales{1}./maxScales1;
        scales{2} = scales{2}./maxScales2;
        x=1:size(scales{1},2); 
        K=0.5; 
        bar1=bar(x, scales{1}, 'FaceColor', 'b', 'EdgeColor', 'b'); 
        set(bar1,'BarWidth',K); 
        hold on;
        bar2=bar(x, scales{2}, 'FaceColor', 'r', 'EdgeColor', 'r');
        set(bar2,'BarWidth',K/2); 
        hold off; 
        legend(['scales1 (x ' num2str(maxScales1) ')'],['scales2 (x ' num2str(maxScales2) ')'])
        scales = origScales;
    else
        for i=1:totalModels
            if totalModels > 1, subplot(1,totalModels, i) ; end
            bar(scales{i}) 
        end
    end
else
    for i=1:totalModels
        if strcmp(model.comp{i}.kern.type, 'rbfardjit')
            fprintf('# Scales of model %d: ', i)
            fprintf('%.4f  ',model.comp{i}.kern.inputScales);
            fprintf('\n');
            scales{i} = model.comp{i}.kern.inputScales;
        else
            fprintf('# Scales of model %d: ', i)
            fprintf('%.4f  ',model.comp{i}.kern.comp{1}.inputScales);
            fprintf('\n');
            scales{i} = model.comp{i}.kern.comp{1}.inputScales;
        end
    end
end
