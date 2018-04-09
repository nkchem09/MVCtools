function [r,rmsecv,yp,rpd_value] = loocv2(x,y,lv,sel)
if nargin <4;sel = 0;end
T = size(x,1);
yp = zeros(size(y));
for i = 1:T
    x1 = x(i,:);
    ind = 1:T;ind(i) = [];
    x2 = x(ind,:);
    y2 = y(ind,:);
    if eq(sel,0)
        [x2,mx] = mncn(x2);
        x1 = scale(x1,mx);
        [y2,my] = mncn(y2);
    elseif eq(sel,1)
        [x2,mx,stdx] = auto(x2);
        x1 = scale(x1,mx,stdx);
        [y2,my,stdy] = auto(y2);
    end
    b = simpls(x2,y2,lv);
    yp0 = x1 * b(:,lv);
    if eq(sel,0)
        yp0 = rescale(yp0,my);
    elseif eq(sel,1)
        yp0 = rescale(yp0,my,stdy);
    end
    yp(i) = yp0;    
end
r = corr(yp,y);
rmsecv = rms(y-yp);
[~,rpd_value] = rpd(yp,y);
function [mcx,mx] = mncn(x)
%MNCN Mean center scales matrix to mean zero.
%  Mean centers matrix (x), returning a matrix with
%  mean zero columns (mcx) and the vector of means
%  (mx) used in the scaling.
%
%I/O: [mcx,mx] = mncn(x);
%
%See also: AUTO, MDAUTO, MDMNCN, MDRESCAL, MDSCALE, SCALE, RESCALE

%Copyright Eigenvector Research 1991-98
%Modified 11/93
%Checked on MATLAB 5 by BMW  1/4/97

[m,n] = size(x);
mx    = mean(x);
mcx   = (x-mx(ones(m,1),:));

function sx = scale(x,means,stds)
%SCALE Scales matrix as specified.
%  Scales a matrix (x) using means (mx) and standard 
%  deviations (stds) specified.
%
%I/O:  sx = scale(x,mx,stdx);
%
%  If only two input arguments are supplied then the function
%  will not do variance scaling, but only vector subtraction.
%
%I/O:  sx = scale(x,mx);
%
%See also: AUTO, MDAUTO, MDMNCN, MDRESCAL, MDSCALE, MNCN, RESCALE

%Copyright Eigenvector Research, Inc. 1991-98
%Modified 11/93 
%Checked on MATLAB 5 by BMW  1/4/97

[m,n] = size(x);
if nargin == 3
  sx = (x-means(ones(m,1),:))./stds(ones(m,1),:);
else
  sx = (x-means(ones(m,1),:));
end
function rx = rescale(x,mx,stdx)
%RESCALE Rescales matrix 
%  Rescales a matrix (x) using the means (mx) and
%  standard deviation (stdx) vectors specified.
%
%I/O: rx = rescale(x,mx,stdx);
%
%  If only two input arguments are supplied the
%  function rescales the means only.
%
%I/O: rx = rescale(x,mx);
%
%See also:  AUTO, MDAUTO, MDMNCN, MDRESCAL, MDSCALE, MNCN, SCALE

%Copyright Eigenvector Research, Inc. 1991-98
%Modified 11/93
%Checked on MATLAB 5 by BMW  1/4/97

[m,n] = size(x);
if nargin == 3
  rx  = (x.*stdx(ones(m,1),:))+mx(ones(m,1),:);
else
  rx  = x+mx(ones(m,1),:);
end

function [ax,mx,stdx] = auto(x)
%AUTO Autoscales matrix to mean zero unit variance
%  Autoscales a matrix (x) and returns the resulting matrix (ax)
%  with mean-zero unit variance columns, a vector of means (mx) 
%  and a vector of standard deviations (stdx) used in the scaling.
%
%I/O:  [ax,mx,stdx] = auto(x);
%
%See also: MDAUTO, MDMNCN, MDRESCAL, MDSCALE, MNCN, SCALE, RESCALE

%Copyright Eigenvector Research, Inc. 1991-98
%Modified 11/93
%Checked on MATLAB 5 by BMW  1/4/97

[m,n] = size(x);
mx    = mean(x);
stdx  = std(x);
ax    = (x-mx(ones(m,1),:))./stdx(ones(m,1),:);

function [SEP,RPD,relSEP] = rpd(xp,xm)

SEP = sqrt(sum((xp-xm).^2) ./ (length(xp)));
relSEP = SEP/mean(xm)*100;
RPD = std(xm) ./ SEP;