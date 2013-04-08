
%-- interpolation
[Yall, lbls] = vargplvmLoadData(dataSetName);
height = lbls(1); width = lbls(2);

timeStampsAll = sort([timeStampsTraining; timeStampsTest]);

dt = (timeStampsAll(end) - timeStampsAll(end-1));
dtt = dt/4;

interp = [timeStampsAll(1):dtt:timeStampsAll(end)]';

timeStampsAll = round(timeStampsAll*100000)/100000;
interp = round(interp*100000)/100000;

timeStampsInterp = setdiff(interp, timeStampsAll);

% The following should give a matrix full of dtt's
% diff(sort([timeStampsAll; timeStampsInterp]))';


[Testmeans2 Testcovars2] = vargplvmPredictPoint(model.dynamics, interp);
Varmu2 = vargplvmPosteriorMeanVar(model, Testmeans2, Testcovars2);

playMov(height, width, 0.001, Varmu2);
% close
% i=1; jOld=1;
% pp = 0.03;
% 
% for i=1:size(Yall,1)
%     imagesc(reshape(Yall(i,:), height, width)); colormap('gray')
%     pause(pp)
%     %fprintf('\n%d\n',i)
%     for j=jOld:jOld+(dt/dtt)-2
%         %fprintf('%d ',j)
%         imagesc(reshape(Varmu2(j,:), height, width));
%         pause(pp)
%     end
%     jOld = j+1;
% end
% 
% 
