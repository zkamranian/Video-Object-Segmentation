% Written by Zahra Kamranian, 2019
% Copyright (c) 2019, Zahra Kamranian, University of Isfahan (zahra.kamranian@eng.ui.ac.ir)
% All rights reserved
% 
% This Function is to create initial foreground object and background regions.
% input:
% 'heatmap_deconv' = image heatmap.
% 'net', here is VGG-16
% 'im' = the input image
% 'im_' is the pre-processed image for VGG16.
% 'seg', is the over-segmentation of the input image. 
% 'flowMap' is the motion map.
% 'labelss' is the oversegmntaion of the image into regions.
% 'img_names' = the name of the image
% 'net', here is VGG-16
% 'featureVector' is the regions features.
% outputs:
% 'scribbled_im' is the image with some initial foreground and background,
% declared by red and green.


function scribble_im = create_foreg_bkg_in_frames_V1(heatmap_deconv,net,im,im_,segs,flowMap,labelss,img_names,featureSeg,param)
       
