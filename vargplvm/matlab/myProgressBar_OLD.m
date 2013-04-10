function pb = myProgressBar_OLD(arg1, arg2, counterPrintStep)

if isstruct(arg1)
    pb = arg1;
    if nargin<2
        return
    end
    if mod(arg2,pb.counterPrintStep) == 0
        fprintf('.');
    end
else
    pb.totalIters = arg1;
    if nargin < 2
        pb.counterLength = 40;
    else
        pb.counterLength = arg2;
    end
    if nargin < 3
        pb.counterPrintStep = round(pb.totalIters/pb.counterLength);
    else
        pb.counterPrintStep = counterPrintStep;
    end
    
    
    for i=1:pb.counterLength
        fprintf(1,'.');
    end
    fprintf(1,'|\n');
end