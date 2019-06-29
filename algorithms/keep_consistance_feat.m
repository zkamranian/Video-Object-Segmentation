% Written by Zahra Kamranian, 2019
% Copyright (c) 2019, Zahra Kamranian, University of Isfahan (zahra.kamranian@eng.ui.ac.ir)
% All rights reserved
%
% This Function is to keep consistance features among all the
% images/frames. and, remove the inconsistent ones.
% inputs:
% 'x' is the feature of all scribbled fg or bkg
% 'hist', is the histograms of all scribbled fg or bg
% outputs:
% 'x' is the refined feature of all scribbled fg or bkg
% 'hsit' is the refined histogram of all scribbled fg or bg
% 'outline' is out features

function [x,hist,outline] = keep_consistance_feat(x,hist)

if nargin < 2, hist = []; end

k = 8;
dist = pdist2(x,x);
for i=1:size(x,1)   
    [sortedX,sortingIndices] = sort(dist(i,:),'descend');
    NC(i,:) = sortingIndices(:,1:k);
end

unique_NC = unique(NC);
out = sortrows([unique_NC,histc(NC(:),unique_NC)],2);
x(out(end-k:end,1),:) = [];

if ~isempty(hist)
    hist(:,out(end-k:end,1)) = [];
end
outline = out(end-k:end,1); 