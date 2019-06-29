% Written by Zahra Kamranian, 2019
% Copyright (c) 2019, Zahra Kamranian, University of Isfahan (zahra.kamranian@eng.ui.ac.ir)
% All rights reserved.
%
% This Function is to extract the hypercolumn features from the regions of
% an input image. The use of GPU accelerates the execution of the function.
% Usually, for the co-segmentation, the hypercolmn features from the fifth 
% layer (and the second layer) obtains best results. However, for other 
% applications (or datasets) you can use the other layers' outputs, as well.
% Inputs: 
% 'im_' is the pre-processed version of an input image. The pre-processing is done
% according to the VGG-16 input.
% 'net' is a convolutional network, here is the VGG-16.
% 'segs' is the over-segmentation of the input image. In the case of
% extracting hypercolumn features for the 'pixel', it would be = [].
% 'pixel_flag' : 0, the hypercolumn features are avaraged for the pixel of the regions.
%             or 1, the hypercolumn features are calculated on th 'pixel' 
% output:
% 'featureVector' is the hypercolumn features for the regions(or pixels) of the input
% image.


function featureVector = extract_hypercolumnF(im_,net,segs,pixel,pixel_flag)

%% feed the input image to the net, and get the results in 'res'
net = vl_simplenn_move(net,'gpu');
res = vl_simplenn(net, im_) ;

%% extracting feature map from each layer and resize it to the input image
%first layer
% FirstConvF1 = gather(res(3).x);
% FirstF1 = imresize(FirstConvF1,net.meta.normalization.imageSize(1:2),'bilinear'); 
% FirstF1 = reshape(FirstF1,[size(FirstF1,1)*size(FirstF1,2),size(FirstF1,3)]);
% FirstF1 = mean(FirstF1,2);
% 
% FirstConvF2 = gather(res(5).x);
% FirstF2 = imresize(FirstConvF2,net.meta.normalization.imageSize(1:2),'bilinear'); 
% FirstF2 = reshape(FirstF2,[size(FirstF2,1)*size(FirstF2,2),size(FirstF2,3)]);
% FirstF2 = mean(FirstF2,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% second layer 
secondConvF1 = gather(res(8).x);
secondF1 = imresize(secondConvF1,net.meta.normalization.imageSize(1:2),'bilinear'); 
secondF1 = reshape(secondF1,[size(secondF1,1)*size(secondF1,2),size(secondF1,3)]);
secondF1 = mean(secondF1,2);

secondConvF2 = gather(res(10).x);
secondF2 = imresize(secondConvF2,net.meta.normalization.imageSize(1:2),'bilinear');
secondF2 = reshape(secondF2,[size(secondF2,1)*size(secondF2,2),size(secondF2,3)]);
secondF2 = mean(secondF2,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% third layer
% thirdConvF1 = gather(res(13).x);
% thirdF1 = imresize(thirdConvF1,net.meta.normalization.imageSize(1:2),'bilinear');
% thirdF1 = reshape(thirdF1,[size(thirdF1,1)*size(thirdF1,2),size(thirdF1,3)]);
% thirdF1 = mean(thirdF1,2);

% thirdConvF2 = gather(res(15).x);
% thirdF2 = imresize(thirdConvF2,net.meta.normalization.imageSize(1:2),'bilinear');
% thirdF2 = reshape(thirdF2,[size(thirdF2,1)*size(thirdF2,2),size(thirdF2,3)]);
% thirdF2 = mean(thirdF2,2);

% thirdConvF3 = gather(res(17).x);
% thirdF3 = imresize(thirdConvF3,net.meta.normalization.imageSize(1:2),'bilinear');
% thirdF3 = reshape(thirdF3,[size(thirdF3,1)*size(thirdF3,2),size(thirdF3,3)]);
% thirdF3 = mean(thirdF3,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fourth layer
% FourthConvF1 = gather(res(20).x);
% FourthF1 = imresize(FourthConvF1,net.meta.normalization.imageSize(1:2),'bilinear');
% FourthF1 = reshape(FourthF1,[size(FourthF1,1)*size(FourthF1,2),size(FourthF1,3)]);

% FourthConvF2 = gather(res(22).x);
% FourthF2 = imresize(FourthConvF2,net.meta.normalization.imageSize(1:2),'bilinear');
% FourthF2 = reshape(FourthF2,[size(FourthF2,1)*size(FourthF2,2),size(FourthF2,3)]);
% FourthF2 = mean(FourthF2,2);

% FourthConvF3 = gather(res(24).x);
% FourthF3 = imresize(FourthConvF3,net.meta.normalization.imageSize(1:2),'bilinear');
% FourthF3 = reshape(FourthF3,[size(FourthF3,1)*size(FourthF3,2),size(FourthF3,3)]);
% FourthF3 = mean(FourthF3,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fifth layer
% FifthConvF1 = gather(res(27).x);
% FifthF1 = imresize(FifthConvF1,net.meta.normalization.imageSize(1:2),'bilinear');
% FifthF1 = reshape(FifthF1,[size(FifthF1,1)*size(FifthF1,2),size(FifthF1,3)]);
% FifthF1 = mean(FifthF1,2);

FifthConvF2 = gather(res(29).x);
FifthF2 = imresize(FifthConvF2,net.meta.normalization.imageSize(1:2),'bilinear');
FifthF2 = reshape(FifthF2,[size(FifthF2,1)*size(FifthF2,2),size(FifthF2,3)]);
FifthF2 = mean(FifthF2,2);

FifthConvF3 = gather(res(31).x);
FifthF3 = imresize(FifthConvF3,net.meta.normalization.imageSize(1:2),'bilinear');
FifthF3 = reshape(FifthF3,[size(FifthF3,1)*size(FifthF3,2),size(FifthF3,3)]);
FifthF3 = mean(FifthF3,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if pixel_flag
   featureVector = double([
                            %squeeze(secondF2(pixel)),... 
                            %squeeze(FifthF2(pixel)),... 
                            squeeze(FifthF3(pixel))...
            ]); 
else
    for s = 1:size(segs,2)
             featureVector(s,:) = double([
                                      %mean(squeeze(secondF2(segs{1,s},:))),... 
                                      %mean(squeeze(FifthF2(segs{1,s},:))),... 
                                      mean(squeeze(FifthF3(segs{1,s},:)))...
                            ]);

    end
end

