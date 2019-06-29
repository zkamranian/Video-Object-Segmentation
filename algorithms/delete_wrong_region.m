% Written by Zahra Kamranian, 2019
% Copyright (c) 2019, Zahra Kamranian, University of Isfahan (zahra.kamranian@eng.ui.ac.ir)
% All rights reserved
% 
% This function is to delete the wrong fg/bkg regions, which are not
% helpfull for co-segmentation.

function [alpha,spLabelInd] = delete_wrong_region(alpha,spLabelInd)

wrongFgInd = [];
wrongBgInd = [];

for ii = 1:size(alpha.fg,1)
    if alpha.fg(ii) < 0.5
        wrongFgInd = [wrongFgInd,ii];        
    end
end
extraFgSp = spLabelInd.fg(wrongFgInd);
extraFgAlpha = alpha.fg(wrongFgInd);
spLabelInd.fg(wrongFgInd) = [];
alpha.fg(wrongFgInd) = [];


for ii = 1:size(alpha.bg,1)
    if alpha.bg(ii) >= 0.5
        wrongBgInd = [wrongBgInd,ii];        
    end
end
extraBgSp = spLabelInd.bg(wrongBgInd);
extraBgAlpha = alpha.bg(wrongBgInd);
spLabelInd.bg(wrongBgInd) = [];
alpha.bg(wrongBgInd) = [];

