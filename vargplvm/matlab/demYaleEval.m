v = 1; % This can be either 1 or 2
modelVis = model.comp{v};
%bar(model.comp{v}.kern.comp{1}.inputScales);
figure
% The following causes OUTOFMEMORY exception except from when we prune the
% video dimensions: (Note: also, lvmVisualise was changed a bit so that no
% dynamic slides is presented, because otherwise a strange error occurs).
modelVis.y = Ytr{v};
sc = 2; % Set to 4 and try again, if you run out of memory
[modelP, newHeight, newWidth] = vargplvmReduceVidModel(modelVis, height, width, sc,sc);
lvmVisualiseGeneral(modelP, [], 'imageMRDVisualise', 'imageMRDModify', false, [newHeight newWidth],0,0,1);
clear modelP
figure,bar(model.comp{1}.kern.comp{1}.inputScales), figure, bar(model.comp{2}.kern.comp{1}.inputScales)

%%
%{
for i=1:size(Ytr{1},1)
curInd = i;
subplot(2,1,1)
         imagesc(reshape(Ytr{1}(curInd,:),height,width)), title(['Ytr1: ' num2str(curInd) ')']), colormap('gray')
subplot(2,1,2)
         imagesc(reshape(Ytr{2}(curInd,:),height,width)), title(['Ytsr: ' num2str(curInd) ')']), colormap('gray')
pause
i
end
%}
%%
%40->57
close
%sharedDims = [1 2 3];
ind1 = 3; %40,57)
ind2 = 54;

x1 = model.X(ind1,:);
x2 = model.X(ind2,:);

y1 = model.comp{1}.y(ind1,:);
y2 = model.comp{1}.y(ind2,:);

z1 = vargplvmPosteriorMeanVar(model.comp{1}, x1);
z2 = vargplvmPosteriorMeanVar(model.comp{1}, x2);

xnew = x1;
xnew(sharedDims) = x2(sharedDims);

z1new = vargplvmPosteriorMeanVar(model.comp{1},xnew);

%%

subplot(2,2,1)
imagesc(reshape(y1,height,width)), title('y1'), colormap('gray')
subplot(2,2,2)
imagesc(reshape(y2,height,width)), title('y2'), colormap('gray')
subplot(2,2,3)
imagesc(reshape(z1,height,width)), title('z1'), colormap('gray')
subplot(2,2,4)
imagesc(reshape(z1new,height,width)), title('z1new'), colormap('gray')

%%
close
for j=1:length(sharedDims)
    ind = sharedDims(j); % one of the shared dims
    for i=-8:0.2:8
        xtemp = x1; xtemp(ind) = x1(ind)+i;
        imagesc(reshape(vargplvmPosteriorMeanVar(model.comp{1},xtemp),height,width)), title(num2str(i)), colormap('gray')
        pause(0.1)
    end
    pause
end