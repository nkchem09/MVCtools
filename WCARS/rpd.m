function [SEP,RPD,relSEP] = rpd(xp,xm)

SEP = sqrt(sum((xp-xm).^2) ./ (length(xp)));
relSEP = SEP/mean(xm)*100;
RPD = std(xm) ./ SEP;