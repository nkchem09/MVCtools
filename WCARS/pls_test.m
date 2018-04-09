function [r,rmsep,y,rpd_value] = pls_test(cal,caltar,calv,caltarv,lv,pre_sel)
    
    if nargin < 6
        pre_sel = 0;
    end
    if eq(pre_sel,0)
        [cal,mx] = mncn(cal);
        calv = scale(calv,mx);
        [caltar,my] = mncn(caltar);
    end
    b = simpls(cal,caltar,lv);
    y = calv* b(:,end);
    if eq(pre_sel,0)
        y = y+my;
    end
    %params
    r = corr(y,caltarv);
    rmsep = rms(y-caltarv);
    [~,rpd_value] = rpd(y,caltarv);
end

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
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
end