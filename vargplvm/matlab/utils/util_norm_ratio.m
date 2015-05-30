function [d,d1] = util_norm_ratio(f,g,displ)
if nargin < 3 || isempty(displ)
    displ = false;
end
 d=norm(f - g)/norm(g + f);
d1=norm(f - g,1)/norm(g + f,1); %%
if displ
    fprintf(1,' Norm1 difference: %d\n Norm2 difference: %d\n',d1,d);
end