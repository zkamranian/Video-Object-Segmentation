% Written by Zahra Kamranian, 2019
% Copyright (c) 2019, Zahra Kamranian, University of Isfahan (zahra.kamranian@eng.ui.ac.ir)
% All rights reserved.
%
% This Function is to create heatmap and scribbled frames using the joint motion
% and appreance maps.
% inputs:
% 'seg', is the over-segmentation of the input image. 
% 'inMap' is the motion map.
% 'labels' is the oversegmntaion of the image into regions.
% 'net', here is VGG-16
% 'featureVector' is the regions features.
% outputs:
% 'Lheatmap_deconv' shows the heatmap.
% 'scribbled_im' shows the image with the intial regions.

function [Lheatmap_deconv scribbled_im]=create_scribbled_frame_V1(seg,inMap,labels,net,featureVector,param)

img_path = param.img_path;
imgstyle = param.imgstyle;
img_names = param.img_names;

scribbled_im=cell(1,size(img_names,1));

for ii=1:size(img_names,1)

    im = imread([img_path ,img_names{ii},'.' imgstyle]); 
    im_ = single(im) ; % note: 0-255 range
    im_ = imresize(im_, net.meta.normalization.imageSize(1:2)) ;
    bsxfun(@minus, im_, net.meta.normalization.averageImage);
    im_ = gpuArray(im_);

    heatmap_deconv = simplenn_deconvolution(net, im_);  
    heatmap_deconv = double(rgb2gray(heatmap_deconv));   
    heatmap_deconv = heatmap_deconv*(1/(max(max(heatmap_deconv))));       
    Lheatmap_deconv{ii}=heatmap_deconv;

    scribbled_im{ii} = create_foreg_bkg_in_frames_V1(Lheatmap_deconv{ii},net,im,im_,seg{ii},inMap{ii},labels{ii},img_names{ii},featureVector{ii},param);  
    
    if rem(ii,50)==0
       close all
    end        
end
