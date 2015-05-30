function fAll = util_artificialShadow(f, step, pad, Sigma, mult, disp)

if nargin < 2 || isempty(step), step=4; end
if nargin < 3 || isempty(pad), pad=8; end
if nargin < 4 || isempty(Sigma), Sigma = [1 0.2; 0.2 0.7]; end
if nargin < 5 || isempty(mult), mult = 2000; end
if nargin < 6 || isempty(disp), disp=false; end

Sigma = Sigma .* mult;

%f = reshape(model.comp{1}.y(1,:), height, width);
%
%f = reshape(model.comp{1}.y(1,200:299),10,10);

% f(f==0)=1;
% mask = f./f;
% 
minf = min(min(f));
maxf = max(max(f));
% ff = (f - minVal) ./ ( maxVal - minVal );
% %your_original_data = minVal + norm_data.*(maxVal - minVal)
% 
% imagesc(mask)

N = size(f,1);
D = size(f,2);

nn = min(N,D);

vv1=1:step:nn;
vv2=nn:-step:1;

fAll = zeros(length(vv1)+length(vv2),size(f,1)*size(f,2));
c=1;
for i=vv1
    yy = sqrt((nn./2)^2 - (i-(nn./2))^2) + nn./2;
    mu = [i yy-pad];
    %Sigma = [.9 .4; .4 .3]'.*1000;
    xx1 = linspace(1,nn,N)';
    xx2 = linspace(1,nn,D)';
    
    [X1,X2] = meshgrid(xx1, xx2);
    X = [X1(:) X2(:)];
    p = mvnpdf(X, mu, Sigma);
    mask = reshape(p,D,N).*10;
    mask = scaleBetween(mask, 0, 255);
    
    %surf(X1,X2,mask);
    
    ff = f+mask';
    fAll(c,:) = ff(:)';
    c = c+1;
    % subplot(1,2,1);
    % imagesc(mask); colormap('gray'); colorbar
    % subplot(1,2,2);
    % imagesc([f ff]); colormap('gray'); colorbar
    if disp
    imagesc(ff); colormap('gray')
    if i==1
        pause
    else
        pause(0.1)
    end
    end
end


for i=vv2
    yy = - sqrt((N./2)^2 - (i-(N./2))^2) + N./2;
    mu = [i yy+pad];
    %Sigma = [.9 .4; .4 .3]'.*1000;
    xx1 = linspace(1,N,N)';
    xx2 = linspace(1,N,D)';
    
    [X1,X2] = meshgrid(xx1, xx2);
    X = [X1(:) X2(:)];
    p = mvnpdf(X, mu, Sigma);
    mask = reshape(p,D,N).*10;
    mask = scaleBetween(mask, 0, 255);
    
    %surf(X1,X2,mask);
    
    ff = f+mask';
    fAll(c,:) = ff(:)';
    c = c+1;
    % subplot(1,2,1);
    % imagesc(mask); colormap('gray'); colorbar
    % subplot(1,2,2);
    % imagesc([f ff]); colormap('gray'); colorbar
    if disp
        imagesc(ff); colormap('gray')
        pause(0.1)
    end
end