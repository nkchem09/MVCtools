clear
%addpath('libPLS_1.95\')
load corn_m51.mat;wavenumber = 1:size(X,2);

spec = X;
targ = y;
% group 
[cal,caltar,calv,caltarv] = ks_group(spec,targ);
lv=10;
pre_method='center';
pre_sel =1;
%none,delbias,dt,snv,msc,1st,1st-snv,snv-1st
param.wid=17;
param.dw = abs(mean(diff(wavenumber)));
param.mspec = mean(cal);
T = 1;
sta_re = cell(T,3);

icars_num_list =100; 
for i = 1:T
    sub_cal = cal;
    sub_caltar = caltar;
    pre_cal = sub_cal;
    pre_calv = calv;

    sel_var = icars_model(pre_cal,sub_caltar,lv,icars_num_list);
    [sub_rcv,sub_rmsecv,~,sub_rpd] = loocv2(pre_cal(:,sel_var),sub_caltar,min(lv,length(sel_var)));
    [sub_r,sub_rmsep,~,sub_rpdv] = pls_test(pre_cal(:,sel_var),sub_caltar,pre_calv(:,sel_var),caltarv,lv);

    sta_re(i,:) = {[sub_rmsecv,sub_rmsep],[sub_rcv,sub_r],sel_var};
    sprintf('RMSECV: %0.6f vs. RMSEP %0.6f',sub_rmsecv,sub_rmsep)
end
%rmpath('libPLS_1.95\')


function v_sel = icars_model(spec,targ,lv,num)
    objs = wcars(num);
    objs.disp_sel = false;
    objs.model(spec,targ,lv);
    v_sel = objs.v_sel;
end

function [cal,caltar,calv,caltarv] = ks_group(spec,targ,pca_sel)
    if nargin < 3
        pca_sel = 'undo';
    end
    if strcmp(pca_sel,'do')
        [~,score] = pca(spec,'NumComponents',12);
    else
        score = spec;
    end
    index=ks(score);
    m = size(score,1);
    m1 = floor(m/3*2);

    model = index(1:m1);
    test = index(m1+1:end);
    cal = spec(model,:);
    calv = spec(test,:);
    caltar = targ(model,:);
    caltarv = targ(test,:);
end
