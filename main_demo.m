% Copyright (c) 2019, Zahra Kamranian, University of Isfahan (zahra.kamranian@eng.ui.ac.ir)
% All rights reserved.


clear all;close all;

addpath 'algorithms'
addpath 'others'
% load the pre-trained CNN
% net = load('imagenet-vgg-verydeep-16.mat') ;

%% *******settings
full_connect = 0; 

% parameters of cosegmentation with local spline regression
param.lambda0 = 1e1;
param.lambda1 = 1e0;
param.lambda2 = 1e-7;%0;%
param.lambda3 = 1e0;
param.gamma = 1e4;
param.ncenter = 1;

dataset = 'bear/';
img_path = ['Datasets/images/Davis/480p/',dataset]; 
Gr_path = ['Datasets/images/Davis/GroundTruth/',dataset];
scribbles_path = ['Datasets/scribbles/Davis/480p/',dataset];
out_path = ['results/Davis/480p/',dataset];

if ~exist(out_path,'file')
    mkdir(out_path);
    mkdir([out_path, 'res']);
end


if ~exist(['Datasets/scribbles/',dataset,'heatmap/'],'file')
    mkdir(['Datasets/scribbles/',dataset, 'heatmap']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read image/frame names
imgstyle = 'jpg' ;
img_dir = dir([img_path '*.' imgstyle]);

n_img = length(img_dir);
img_names_t = cell(n_img,1);
for i =1:n_img
    img_names_t{i} = strtok( img_dir(i).name,'.');  
end

param.img_path = img_path;
param.scribbles_path = scribbles_path;
param.out_path = out_path;
param.Gr_path = Gr_path;
param.img_names = img_names_t;
param.imgstyle = imgstyle;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% *******compute flow and inMaps*******
fprintf('Compute optical flow....\n');
video_path = [img_path,'newVideo/'];
if ~exist(video_path,'file')
    mkdir(video_path);
end

options.infolder = video_path;
options.outfolder = out_path;
options.flowmethod = 'broxPAMI2011';
options.visualise = true;
options.vocal = true;
options.ranges = [ 1, n_img+1 ];
options.positiveRanges = [ 1, 2 ];
options.maxedge = 400;

 for i = 1:n_img
        im = imread([img_path ,img_names_t{i},'.' imgstyle]); 
        im = imresize(im,net.meta.normalization.imageSize(1:2),'nearest');
        imwrite(im,[video_path,img_names_t{i},'.',imgstyle]);
 end

flow = loadFlow( options, 1 );
if isempty( flow ) 
    flow = computeOpticalFlow( options, 1 );
end
inMaps = getInOutMaps( flow );

inMaps{size(img_names_t,1)} = inMaps{size(img_names_t,1)-1};
for i = 1: n_img
    if i>1 & i<n_img
        inMap{i} = inMaps{i-1}+inMaps{i};
    end   
end
inMap{1} = inMap{2};
inMap{n_img} = inMap{n_img-1};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% *******oversegmentation based on ucm algorithm (coarse segmentation)
fprintf('Coarse segmentation of the images/frames....\n');
coarseseg_data05 = [out_path 'coarseseg_data_new05.mat'];
if exist(coarseseg_data05,'file') % & exist(coarseseg_data,'file')
    load(coarseseg_data05);
else
    for i = 1:n_img
        im = imread([img_path ,img_names_t{i},'.' imgstyle]);  
        [coarseLabel05{i} coarseSeg05{i} histSP05{i} colors05{i} edges05{i}]=coarse_Seg(im,full_connect,thr);      
        if rem(i,50)==0
            close all
        end
    end
    save(coarseseg_data05,'coarseLabel05','coarseSeg05','histSP05', 'colors05','edges05');  
end
close all

%% *******extract hypercolumn features of the regions
fprintf('Extracting hypercolumn features from the regions....\n');
Feature_VecFile = [out_path 'FeatureVec.mat'];
if exist(Feature_VecFile,'file') 
     load(Feature_VecFile);
else
    for ii=1:size(img_names_t,1)
        im = imread([img_path ,img_names_t{ii},'.' imgstyle]); 
        im_ = single(im) ; % note: 0-255 range
        im_ = imresize(im_, net.meta.normalization.imageSize(1:2)) ;
        bsxfun(@minus, im_, net.meta.normalization.averageImage);
        im_ = gpuArray(im_);
        featureVector{ii} = extract_hypercolumnF(im_,net,coarseSeg05{ii},[],0);
    end
    save(Feature_VecFile,'featureVector');
end

%% *******create scribbled images/frames
fprintf('Creating heatmap from the imgages/frames....\n');
Heat_mapFile=[out_path 'heatmapSP_coarse_new05_sp.mat'];
if exist(Heat_mapFile,'file') 
    load(Heat_mapFile)
else
      % create scribbled frames
    [Lheatmap, scribbled_im] = create_scribbled_frame_V1(coarseSeg05,inMap,coarseLabel05,net,featureVector,param);
    save(Heat_mapFile, 'Lheatmap','scribbled_im');
end
close all;

%% *******read the scribbled image names
scribbles_dir = dir([scribbles_path '*.bmp']);
n_scri = length(scribbles_dir);
scribbles_names = cell(n_scri,1);
scri_img_idx = zeros(n_scri,1);%index of scribbles image
for i =1:n_scri
    temp_name =  strtok(scribbles_dir(i).name,'.');  
    scribbles_names{i} = temp_name;
    j=1;
    for j = 1:n_img
        if strcmp(temp_name,img_names_t{j})
            scri_img_idx(i) = j;
            break
        end
    end
    
end
%% *******resort image name ,scribbled image in the front
un_scri_img_idx = ones(n_img,1);
un_scri_img_idx(scri_img_idx) = 0;
un_scri_img_idx = find(un_scri_img_idx);
re_img_idx = [scri_img_idx;un_scri_img_idx];
img_names = cell(n_img,1); histSP = cell(n_img,1);
labels = cell(n_img,1); colors_s = cell(n_img,1);
lab_colors_s = cell(n_img,1); edges_s = cell(n_img,1);
seg = cell(n_img,1); d_edges = cell(n_img,1);
for i  = 1:n_img
      img_names{i} = img_names_t{re_img_idx(i)};
      histSP{i} = histSP05{re_img_idx(i)};
      labels{i} = coarseLabel05{re_img_idx(i)};
      colors_s{i} = colors05{re_img_idx(i)};
      edges_s{i} = edges05{re_img_idx(i)};
      seg{i} = coarseSeg05{re_img_idx(i)};
end
clear un_scri_img_idx  re_img_idx 

%% *******get scribbles label of superpixel
superpixel_labelInd = cell(n_scri,1);
for i = 1:n_scri
    scribs_img_name = [scribbles_path scribbles_names{i} '.bmp'];
    [lines] = seed_generation_1(scribs_img_name);
    fg = unique(labels{i}(find(lines(:,1))));
    bg = unique(labels{i}(find(lines(:,2))));
    tmp1 = [fg;bg]; nf = length(fg);
    [b1 m1 n1] = unique(tmp1,'first');
    temp_label.fg = b1(m1<=nf); temp_label.bg = b1(m1>nf);
    superpixel_labelInd{i} = temp_label;
end

%% *******Co-segmentation
fprintf('Cosegmentation...\n');

[YS,res] = Co_segmentation_V1(histSP,superpixel_labelInd,edges_s,featureVector,labels,seg,Lheatmap,param);

% transform the superpixels results to pixels results 
for i = 1:n_img
    YS_temp = YS{i};
    [h,w,d] = size(labels{i});
    ind = find(YS_temp);
    Y = zeros(h,w);
    for j = 1:length(ind)
        Y(labels{i} == ind(j))=1;
    end
    Y_mask = uint8(repmat(Y,[1,1,3]));
    im = imread([img_path ,img_names{i},'.' imgstyle]);
    im = imresize(im,[224 224],'nearest');

   [imgMasks,segOutline,imgMarkup]=segoutput(im2double(im),double(Y+1));

    % save pixels results
    file_save = [out_path, 'res/',img_names{i}, '_segmentation.jpg'];
    file_save_img = [out_path, 'res/',img_names{i}, '_mask_img.jpg'];
    imwrite(Y,  file_save);
    imwrite(imgMarkup,  file_save_img);
    unique(Y);
       
    
    %% iCoseg
%     groundtruth_path = [img_path,'GroundTruth/'];
%     gtImage = imread([groundtruth_path,img_names_t{i},'.png']);
%     groundtruth = double(gtImage(:,:,1)>0); groundtruth=imresize(groundtruth,[224 224],'nearest'); %unique(groundtruth); imshow(groundtruth);
%      %imshow(groundtruth);
%     P(i) =sum(groundtruth(:)==Y(:)) ./ prod(size(groundtruth));
%     Jar(i) =sum( (Y(:)==1) & (groundtruth(:)==1) ) ./ sum( (Y(:) | groundtruth(:))==1 );

    %% MSRC/Davis/SegTrackV2
    gtImage = imread([Gr_path,img_names{i},'.png']);  
    gtImage(gtImage>200)=255; gtImage(gtImage<=200)=0;      
    groundtruth = double(gtImage(:,:,1))./255; unique(groundtruth);
    groundtruth=imresize(groundtruth,[224 224],'nearest'); %unique(groundtruth); imshow(groundtruth);
    P(i) =sum(groundtruth(:)==Y(:)) ./ prod(size(groundtruth));
    if sum((groundtruth(:)==1) )
        Jar(i) =sum( (Y(:)==1) & (groundtruth(:)==1) ) ./ sum( (Y(:) | groundtruth(:))==1 );
    else
        if sum((Y(:)==1))
            Jar(i) = 0;
        else
            Jar(i) =1;
        end            
    end
  %% Internet Images
%     groundtruth_path = Gr_path;
%     gtImage = imread([groundtruth_path,img_names{i},'.png']);   
%     groundtruth = double(gtImage(:,:,1)>0); groundtruth=imresize(groundtruth,[224 224],'nearest'); %unique(groundtruth); imshow(groundtruth);
%     P(i) =sum(groundtruth(:)==Y(:)) ./ prod(size(groundtruth));
%     Jar(i) =sum( (Y(:)==1) & (groundtruth(:)==1) ) ./ sum( (Y(:) | groundtruth(:))==1 );
%     Precision(i) =sum( (Y(:)==1) & (groundtruth(:)==1) ) ./ (sum( (Y(:)==1) & (groundtruth(:)==1) )+sum( (Y(:)==1) & (groundtruth(:)==0) ));
end


MP = mean(P); MJ = mean(Jar); 
save ([out_path,'res/PJ'] ,'P','Jar','MP','MJ','img_names');
fprintf('Accuracy = %f\n   IOU   = %f\n',MP,MJ);



