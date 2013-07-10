function h = vargplvmPrintPlot(model, lbls, capName, experimentNo, writeFile, skipPng)

% VARGPLVMPRINTPLOT Print latent space for learnt model.
%
% h = vargplvmPrintPlot(model, lbls, capName, experimentNo)
% VARGPLVM

if nargin < 6
    skipPng = false;
end

if nargin < 5 || isempty(writeFile)
    writeFile = false;
end


lvmScatterPlot(model, lbls);

if writeFile
    fileName = ['dem' capName num2str(experimentNo)];
    fprintf('# Printing %s...\n', fileName);
    try
        print('-depsc', ['../tex/diagrams/' fileName])
        print('-deps', ['../tex/diagrams/' fileName 'NoColour'])
    catch
        mkdir('diagrams');
        print('-depsc', ['./diagrams/' fileName])
        print('-deps', ['./diagrams/' fileName 'NoColour'])
    end
end

if skipPng
    return
end




% make smaller for PNG plot.
pos = get(gcf, 'paperposition');
origpos = pos;
pos(3) = pos(3)/2;
pos(4) = pos(4)/2;
set(gcf, 'paperposition', pos);
fontsize = get(gca, 'fontsize');
set(gca, 'fontsize', fontsize/2);
lineWidth = get(gca, 'lineWidth');
set(gca, 'lineWidth', lineWidth*2);
if writeFile
    try
        print('-dpng', ['../html/' fileName])
    catch
        mkdir('html');
        print('-dpng', ['./html/' fileName])
    end
end
set(gcf, 'paperposition', origpos);

figure
clf
ax = axes('position', [0.05 0.05 0.9 0.9]);
hold on
lvmTwoDPlot(model.X, lbls, getSymbols(size(lbls, 2)));
xLim = [min(model.X(:, 1)) max(model.X(:, 1))]*1.1;
yLim = [min(model.X(:, 2)) max(model.X(:, 2))]*1.1;
set(ax, 'xLim', xLim);
set(ax, 'yLim', yLim);
set(ax, 'fontname', 'arial');
set(ax, 'fontsize', 20);
if writeFile
    try
        print('-depsc', ['../tex/diagrams/' fileName 'NoGray'])
        print('-deps', ['../tex/diagrams/' fileName 'NoGrayNoColour'])
    catch
        print('-depsc', ['./diagrams/' fileName 'NoGray'])
        print('-deps', ['./diagrams/' fileName 'NoGrayNoColour'])
    end
end
