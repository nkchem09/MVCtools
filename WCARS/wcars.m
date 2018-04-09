%{
update max(selecting run) set 100
一种基于窗口竞争性自适应重加权采样策略的近红外特征变量选择方法
code:Guorong Du, Pao Li
%}
classdef wcars < handle
    properties
        fold_num
        i_num     
        mc_ratio = 0.8;
        pre_method = 'center';%autoscaling
        b_r
        para_x_r
        para_y_r
        b_all
        w_all
        v_sel
        cut_num
        rmsecv_r
        disp_sel = true;
    end
    methods
        %fun begin
        %define
        function obj = wcars(i_num,fold_num,pre_method)
            if nargin > 0
                obj.i_num = i_num;
            else
                obj.i_num = 100;
            end
            if nargin > 1
                obj.fold_num = fold_num;
            else
                obj.fold_num = 10;
            end            
            if nargin > 2
                obj.pre_method = pre_method;
            else
                obj.pre_method = 'center';
            end
        end
        %fun begin
        function obj = model(obj,spec,targ,lv)
            spec0 = spec;
            [m,n] = size(spec);
            vs    = floor(m * obj.mc_ratio);
            wnum = obj.i_num;
            N = min(100,wnum-1);
            ratio1 = (wnum-1)/N;
            
            wid  = floor(linspace(1,n,wnum+1));
            wid2 = [wid(1:end-1)',wid(2:end)'-1];
            wid2(1,1) = 1;
            wid2(end,2) = n;
            
            ind2 = ones(n,N);
            B2   = zeros(size(ind2));
            w2   = zeros(size(ind2));
            rmsecv = zeros(N,1);
            
            for i = 1:N
                if obj.disp_sel
                    if eq(i/100,round(i/100))
                        sprintf('modeling step: %d out of %d is running',i,N)
                    end
                end
                rmsecv(i)    = obj.kfcross(spec,targ,lv,10); 
                [cal,caltar] = obj.rgroup(spec,targ,vs);
                cal = pre_data_model(cal,obj.pre_method);
                caltar = pre_data_model(caltar,obj.pre_method);
                b = simpls(cal,caltar,min([lv,size(cal,2)]));
                B2(:,i) = b(:,end);
                ab           = abs(b(:,end));
                w2(:,i)      = ab ./sum(ab);
                %  
                l = zeros(wnum,1);
                for j = 1:wnum
                    l(j) = mean(ab(wid2(j,1):wid2(j,2)));
                end
                [~,temp] = sort(l);
                for j = 1:floor(i*ratio1)
                    ind2(wid2(temp(j),1):wid2(temp(j),2),i+1) = 0;
                    spec(:,wid2(temp(j),1):wid2(temp(j),2)) = 0;
                end
            end
            [~,vl] = min(rmsecv);
            vind = eq(ind2(:,vl),1);
            [xc,para_xc] = pre_data_model(spec0(:,vind),obj.pre_method);
            [yc,para_yc] = pre_data_model(targ,obj.pre_method);
            obj.b_r = simpls(xc,yc,min([lv,size(xc,2)]));
            obj.para_x_r = para_xc;
            obj.para_y_r = para_yc;
            obj.b_all = B2;
            obj.w_all = w2;
            obj.v_sel = vind;
            obj.cut_num= vl;
            obj.rmsecv_r = rmsecv;
        end
        %fun begin
        function [yp,r,rmsep] = pred(obj,specv,targv)
            xv = pre_data_pre(specv(:,obj.v_sel),obj.para_x_r,obj.pre_method);
            xv(isinf(xv)) = 0;
            yp = xv*obj.b_r(:,end);
            %rescaling
            yp = pre_data_recovery(yp,obj.para_y_r,obj.pre_method);
            r = corr(yp,targv);
            rmsep = rms(yp-targv);
        end
    end
    %method static
    methods (Static)
        %fun new
        function [cal,caltar,calv,caltarv] = rgroup(spec,targ,n)
        %
        ind       = randperm(length(targ));
        indt      = ind(1:n);
        indv      = ind;
        indv(1:n) = [];
        %
        cal       = spec(indt,:);
        calv      = spec(indv,:);
        caltar    = targ(indt);
        caltarv   = targ(indv);
        end
        %fun new
        function rmsecv = kfcross(spec,targ,lv,k,pre_method)
            if nargin < 4
                k= 10;
            end
            if nargin < 5
                pre_method = 'center';
            end
            m = size(spec,1);
            yp = zeros(size(targ));
            for i = 1:k
                ind = i:k:m;
                calv = spec(ind,:);
                cal = spec;
                cal(ind,:) = [];
                caltar = targ;
                caltar(ind) = [];
                
                if strcmpi(pre_method,'center')
                    [cal,mx] = mncn(cal);
                    calv = scale(calv,mx);
                    [caltar,my]=mncn(caltar);
                elseif strcmpi(pre_method,'autoscaling')
                    [cal,mx,stdx] = zscore(cal);
                    calv = scale(calv,mx,stdx);
                    [caltar,my,stdy] = zscore(caltar);
                else
                    sprintf('Warning: kf-cross validation without pre-processing!')
                end
                b   =  simpls(cal,caltar,lv);
                y   =  calv * b(:,lv);
                if strcmpi(pre_method,'center')
                    y = y + my;
                elseif strcmpi(pre_method,'autoscaling')
                    y = y*stdy+my;
                else
                    sprintf('Warning: kf-cross validation without pre-processing!')
                end
                yp(ind) = y;
            end
            rmsecv = rms(yp - targ);
        end
        %
    end
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

function [x,para] = pre_data_model(x,method)
    if strcmpi(method,'center')
        [x,mx] = mncn(x);
        para.mx = mx;
    elseif strcmpi(method,'autoscaling')
        [x,mx,stdx] = zscore(x);
        para.mx = mx;
        para.stdx = stdx;
    else
        x = x;
        para = [];
    end
end
function x = pre_data_pre(x,para,method)
    if strcmpi(method,'center')
        x = scale(x,para.mx);
    elseif strcmpi(method,'autoscaling')
        x = scale(x,para.mx,para.stdx);
    else
        x = x;
    end
end
function x = pre_data_recovery(x,para,method)
    if strcmpi(method,'center')
        x = rescale(x,para.mx);
    elseif strcmpi(method,'autoscaling')
        x = rescale(x,para.mx,para.stdx);
    else
        x = x;
    end
end
function r = rms(y)
% rms -- Return Root Mean Square of Signal
%  Usage
%    r = rms(y)
%  Inputs
%    y    signal
%  Outputs
%    r    sqrt( mean ( y .^2 ))
%
r = sqrt( mean ( y .^2 ));
%   
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%   
end    
%
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
end

% ------------------------------------------------------------------------
% Function: [B,C,P,T,U,R,R2X,R2Y]=plssim(X,Y,h,S,XtX)
% ------------------------------------------------------------------------
% Aim:
% Partial Least Squares for tall X matrices, SIM-PLS 
% ------------------------------------------------------------------------
% Input: 
% X, matrix (n,p), predictor matrix (assumed to be centered)
% Y, matrix (n,m), predictand (assumed to be centered)
% h, scalar, number of PLS factors
% S, matrix (n,m), S=X'*Y
% XtX, matrix (n,n), XtX=X'*X (boosts speed for tall X matrices n>>p)
% ------------------------------------------------------------------------
% Output: 
% B, matrix (p,m), regression coefficients
% C, matrix (m,h), Y loadings
% P, matrix (p,h), X loadings
% T, matrix (n,h), X scores (standardized) 
% U, matrix (n,h), Y scores
% R, matrix (p,h), X weights
% R2X, vecor (1,h), X-variance
% R2Y, vecor (1,h), Y-variance
% ------------------------------------------------------------------------
% Example: 
% 1/ for non tall matrix: 
% [B]=plssim(X,Y,10,[],X'*Y)
% 2/ for tall matrix:     
% [B]=plssim(X,Y,10,X'*Y,X'*X)
% ------------------------------------------------------------------------
% The above routine is included into the toolbox with personal agreement 
% of its author Sijmen de Jong
% ------------------------------------------------------------------------
% Reference:
% S. de Jong, SIMPLS: An alternative approach to partial least squares 
% regression, Chemometrics and Intelligent Laboratory Systems, 
% 18 (1993) 251-263

function [B,C,P,T,U,R,R2X,R2Y]=simpls(X,Y,A,S,XtX);

[n,px]=size(X); 
[n,m]=size(Y);   				

if nargin<5, 
    S=[]; 
end

if isempty(S) 
    S=(Y'*X)'; 
end		            % if XtX not inputted, S=[]; always when S=[] then S=(Y'*X)'

if nargin<4 
    XtX=[]; 
end					% if S is not inputted, XtX=[];

if isempty(XtX) & n>3*px 
    XtX=X'*X; 
end			        % when XtX=[] and X is very "tall", the booster XtX is calculated

if nargin<3 
    A=10; 
end

A=min([A px n-1]);			% if A is not inputted, then the defaul A is min[10 px n-1]
T=zeros(n,A); 
U=T;						% initialization of variables
R=zeros(px,A); 
P=R; 
V=R;
C=zeros(m,A); 
R2Y=zeros(1,A);
z=zeros(m,1); 
v=zeros(px,1);

if n>px 
    S0 = S; 
end
StS=S'*S;				    % SIMPLS algorithm
nm1=n-1;
tol=0;

for a=1:A
    StS=StS-z*z'; 
    [Q,LAMBDA]=eig(StS); 
    [lambda,j]=max(diag(LAMBDA)); 
    q=Q(:,j(1));
    r=S*q;
    t=X*r;
    
    if isempty(XtX)
        p=(t'*X)'; 
    else
        p=XtX*r;
    end
    
    if n>px, 
        d=sqrt(r'*p/nm1); 
    else 
        d=sqrt(t'*t/nm1); 
    end
    
    v=p-V(:,1:max(1,a-1))*(p'*V(:,1:max(1,a-1)))'; 
    v=v/sqrt(v'*v); 
    z=(v'*S)'; 
    S=S-v*z'; 
    V(:,a)=v;
    R(:,a)=r/d; 						    % X weights
    P(:,a)=p/(d*nm1); 						% X loadings
    T(:,a)=t/d;							    % X scores
    U(:,a)=Y*q;							    % Y scores
    C(:,a)=q*(lambda(1)/(nm1*d)); 			% Y loadings
    R2Y(1,a)=lambda(1)/d;					% Y-variance accounted for
    B(:,a)=R*C';					        % B-coefficients of the regression Y on X
   
end

clear StS V LAMBDA Q p q r t v z;

if d<tol,
    A=a-1; 
    a=A; 
    T=T(:,1:A); 
    U=U(:,1:A); 
    R=R(:,1:A); 
    P=P(:,1:A); 
    C=C(:,1:A);
end

while a>1
    U(:,a) = U(:,a)-T(:,1:a-1)*(U(:,a)'*T(:,1:a-1)/nm1)'; 
    a=a-1; 
end

if isempty(XtX),
    sumX2=sum(X.^2);
else 
    sumX2=sum(diag(XtX)); 
end

R2X=100*nm1/sum(sumX2)*(sum(P.^2)); 
R2Y=100/nm1/sum(sum(Y.^2))*(R2Y(1:A).^2);
end