% Plot points with different color/symbol according to targets t

function pl = plot_cl(y, dims, t, symb)

if nargin < 3 || isempty(t)
    t = ones(size(y,1), 1);
    un_t = 1;
else
    un_t = unique(t);
end

if nargin < 4 || isempty(symb)
    symb = getSymbols(length(un_t));
else
    if length(un_t) == 1 && ~iscell(symb)
        symb = {symb};
    end
end
if length(un_t) > length(symb)
    symb = getSymbols(length(un_t));
end

legnd = {};
pl = {};
for i=1:length(un_t)
    ind = t == un_t(i);
    switch length(dims)
        case 1
            pl{end+1} = plot(y(ind,dims(1)), symb{i});
        case 2
            pl{end+1} = plot(y(ind,dims(1)), y(ind, dims(2)), symb{i});
        case 3
            pl{end+1} = plot3(y(ind, dims(1)), y(ind, dims(2)), y(ind, dims(3)), symb{i});
             grid on
    end
    legnd{end+1} = num2str(un_t(i));
    hold on;
end
legend(legnd)
