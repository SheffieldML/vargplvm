subsets = 100:100:1000;

 svargplvmErrorsLog = [20,18,8,6,3,5,6,6,5,5];
 NNerrorsLog = [42,24,12,9,9,10,10,9,8,8];

 
svargplvmErrors = [];
NNerrors = [];
varianceAll = [];

allExperiments = [3334 3336:3344];
for i=1:length(allExperiments)
    try
    curExp = allExperiments(i);
   % load(['demOilSvargplvm' num2str(curExp)])
   % model.N
    load(['demOilSvargplvmPred' num2str(curExp)])
     if errors.gplvm~=svargplvmErrorsLog(i)
         warning(['Mismatch for svargplvm in experiment ' num2str(curExp)])
     end
      if  errors.NN ~= NNerrorsLog(i)
         warning(['Mismatch for NN in experiment ' num2str(curExp)])
      end
   % fprintf('Exp: %d. ErrorGplvm:%.3f - ErrorGplvmLogs: %.3f\n',curExp, errors.gplvm, svargplvmErrors(i));
   svargplvmErrors = [svargplvmErrors errors.gplvm];
   NNerrors = [NNerrors errors.NN];
     curVar = sum(sum(predictions.varX))/length(predictions.varX(:));
     varianceAll = [varianceAll curVar];
%     catch %%%%%%TEMP
%            svargplvmErrors = [svargplvmErrors 3];
%          NNerrors = [NNerrors 9];
%          curVar = 0;
%         varianceAll = [varianceAll curVar]; %%%%%%%%%TEMP
    end %%%%TEMP
end


svargplvmErrors = svargplvmErrorsLog;%%%%%%%%%
NNerrors = NNerrorsLog;%%%

accuracySvargplvm = (1000 - svargplvmErrors)/10;
accuracyNN = (1000 - NNerrors)/10;

% plot(subsets, svargplvmErrors, 'x-')
% hold on
% plot(subsets, NNerrors, 'x-r')
% hold off

%errorbar(subsets,svargplvmErrors/100,varianceAll/2)

plot(subsets, accuracySvargplvm, '.--', 'MarkerSize', 18)
hold on
plot(subsets, accuracyNN, '.--r','MarkerSize', 18)
hold off
set(gca,'FontSize',12)
legend('MRD','NN','FontSize',12)
xlabel('Number of training points','FontSize',12)
ylabel('Accuracy %','FontSize',12)