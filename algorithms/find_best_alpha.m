% Written by Zahra Kamranian, 2019
% Copyright (c) 2019, Zahra Kamranian, University of Isfahan (zahra.kamranian@eng.ui.ac.ir)
% All rights reserved
% 
% This function is to find the best alpha vector which minimize the cost
% function. During the function, the initial fg/bkg regions (which are
% obtained from the visualization step, but is not helpful for
% co-segmentation) are removed.
% input:
% 'superpixel_labelInd' is the labeled fg/bkg regions
% 'YS' is the segmentation results
% 'featureVec' is the feature vector
% outputs:
% 'superpixel_labelInd' is the refined fg/bkg regions
function [superpixel_labelInd] = find_best_alpha(superpixel_labelInd,YS,featureVec,param)

options.type_interaction = 1; 
options.gamma_global     = param.lambda0;
options.lambda_local     = 0.0001;

n_img = length(YS);
n_scri = length(superpixel_labelInd);


%% scribbles images are the supervised images
n_sp = zeros(n_img,1);%the number of superpixel in each image
superpixel_labelInd_all.fg = [];
superpixel_labelInd_all.bg = [];
M = cell(n_img,1);
S=0;

%get all scribbled superpixels
scribbled_sp_colors_all_fg = [];
scribbled_sp_colors_all_bg = [];
for i = 1:n_img
    scribbled_sp_colors_all_fg = [scribbled_sp_colors_all_fg;featureVec{i}(find(YS{i}==1),:)];
    scribbled_sp_colors_all_bg = [scribbled_sp_colors_all_bg;featureVec{i}(find(YS{i}==0),:)];
end

% remove unconsistance features
[scribbled_sp_colors_all_fg] = keep_consistance_feat(scribbled_sp_colors_all_fg);
[scribbled_sp_colors_all_bg] = keep_consistance_feat(scribbled_sp_colors_all_bg);

scribbled_sp_colors_all = [scribbled_sp_colors_all_fg;scribbled_sp_colors_all_bg]';%d x n
scribbled_sp_all_ind.fg = [1:size(scribbled_sp_colors_all_fg,1)]';
scribbled_sp_all_ind.bg = [size(scribbled_sp_colors_all_fg,1)+1:size(scribbled_sp_colors_all,2)]';

%makeGMM features
A = scribbled_sp_colors_all;
bgPix = scribbled_sp_all_ind.bg ;
fgPix = scribbled_sp_all_ind.fg ;
ncenters = param.ncenter;
GMMparam.ncenters = ncenters;
X2 = A(:,bgPix);
[Unaries, Mu, Sigma, ncentersB] = EM_init_kmeans(X2, ncenters);
[GMMparam.unarypotB, GMMparam.muB, GMMparam.sigmaB] = EM(X2, Unaries, Mu, Sigma);
GMMparam.ncentersB = ncentersB;
X2 = X2';


X1 = A(:,fgPix); 
[Unaries, Mu, Sigma, ncentersF] = EM_init_kmeans(X1, ncenters);
[GMMparam.unarypotF, GMMparam.muF, GMMparam.sigmaF] = EM(X1, Unaries, Mu, Sigma);
GMMparam.ncentersF = ncentersF;
labelParam = GMMparam;
labelParam.fea_style = 'features';
X1 = X1';


for i = 1:n_scri
    initialAlphaFg = ones(size(superpixel_labelInd{i}.fg,1),1);
    initialAlphaBg = ones(size(superpixel_labelInd{i}.bg,1),1);
    alpha1{i}.fg = [];   alpha1{i}.bg = [];
      
    if ~isempty(superpixel_labelInd{i}.fg)
         Y = featureVec{i}(superpixel_labelInd{i}.fg,:);
        for ii = 1:size(Y,1)
            alph(ii) = fminunc(@(alpha)myfunalpha(alpha,Y(ii,:),X1),initialAlphaFg(ii));
        end
        alpha1{i}.fg = 1-abs(1-alph);    %1./alph;
        clear alph
    end
       
    
    if ~isempty(superpixel_labelInd{i}.bg)
        Y = featureVec{i}(superpixel_labelInd{i}.bg,:);
        for ii = 1:size(Y,1)
            alph(ii) = fminunc(@(alpha)myfunalpha(alpha,Y(ii,:),X2),initialAlphaBg(ii));
        end
        alpha1{i}.bg = abs(1-alph);     %(1-(1./alph));
    end
   
    
    [alpha1{i},superpixel_labelInd{i}] = delete_wrong_region(alpha1{i},superpixel_labelInd{i});
end

