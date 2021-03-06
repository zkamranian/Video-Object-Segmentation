%
% This Function is for unsupervised images/frames co-segmentation. 
% This is done base on the scribbled images/frames which are obtained from previous steps.
% inputs:
% 'histSP' is the set of histograms
% 'superpixel_labelInd', superpixel_labelInd{i} shows the fg/bkg regions in
% i-th image
% 'edges_s' is a set of edges
% 'featureVector' is a set of features
% 'labels' is a set of labels
% 'seg' is a set of segments
% 'Lheatmap' is a set of heatmaps
% outputs:
% 'YS' is a set of segmentation results
% 'res' is accuracy and IOU results

function [YS,res] = Co_segmentation_V1(histSP,superpixel_labelInd,edges_s,featureVector,labels,seg,Lheatmap,param)
  
lambda0 = param.lambda0;
lambda1 = param.lambda1;
lambda2 = param.lambda2;
lambda3 = param.lambda3;     
img_path = param.img_path;
img_names = param.img_names;
imgstyle = param.imgstyle;
Gr_path = param.Gr_path;

n_img = length(histSP);
n_scri = length(superpixel_labelInd);

superpixel_labelIndNew = cell(n_img,1);

%% calculate heatmap for superpixels
for i=1:n_img
    heat=Lheatmap{i};
    for jj=1:size(seg{i},2)
        he{i}(1,jj)=mean(heat(seg{i}{jj}));
    end
end


counter = 0;
AA = 1;
BB = 1;


while ~isempty(AA) && ~isempty(BB)
    counter = counter+1;
    fprintf('                                                    count = %d\n',counter);
    %% find the consistense features for scribbles images
    n_sp = zeros(n_img,1);%the number of superpixel in each image
    superpixel_labelInd_all.fg = [];
    superpixel_labelInd_all.bg = [];
    M = cell(n_img,1);
    empty_labelInd.hg=[1];
    empty_labelInd.bg=[1];
    S=0;

    %get all scribbled superpixels
    scribbled_sp_feat_all_fg = [];
    scribbled_sp_feat_all_bg = [];
    scribbled_sp_hist_all_fg = [];
    scribbled_sp_hist_all_bg = [];
    sp_fg = [];
    sp_bg = [];
    for i = 1:n_scri
        scribbled_sp_feat_all_fg = [scribbled_sp_feat_all_fg;featureVector{i}(superpixel_labelInd{i}.fg,:)];
        if ~isempty(superpixel_labelInd{i}.fg)
            sp_fg = [sp_fg;ones(size(superpixel_labelInd{i}.fg,1),1)*i,superpixel_labelInd{i}.fg];
        end
        scribbled_sp_feat_all_bg = [scribbled_sp_feat_all_bg;featureVector{i}(superpixel_labelInd{i}.bg,:)];
        if ~isempty(superpixel_labelInd{i}.bg)
            sp_bg = [sp_bg;ones(size(superpixel_labelInd{i}.bg,1),1)*i,superpixel_labelInd{i}.bg];
        end
        scribbled_sp_hist_all_fg = [scribbled_sp_hist_all_fg,histSP{i}(:,superpixel_labelInd{i}.fg)];
        scribbled_sp_hist_all_bg = [scribbled_sp_hist_all_bg,histSP{i}(:,superpixel_labelInd{i}.bg)];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % remove unconsistence features
    [scribbled_sp_feat_all_fg,scribbled_sp_hist_all_fg,outline_fg] = keep_consistance_feat(scribbled_sp_feat_all_fg,scribbled_sp_hist_all_fg);
    [scribbled_sp_feat_all_bg,scribbled_sp_hist_all_bg,outline_bg] = keep_consistance_feat(scribbled_sp_feat_all_bg,scribbled_sp_hist_all_bg);
    %%%%%%%%
   
    scribbled_sp_hist_all = [scribbled_sp_hist_all_fg,scribbled_sp_hist_all_bg];%d x n
    scribbled_sp_feat_all = [scribbled_sp_feat_all_fg;scribbled_sp_feat_all_bg]';%d x n
    scribbled_sp_all_ind.fg = [1:size(scribbled_sp_feat_all_fg,1)]';
    scribbled_sp_all_ind.bg = [size(scribbled_sp_feat_all_fg,1)+1:size(scribbled_sp_feat_all,2)]';


    %makeGMM features
    A = scribbled_sp_feat_all;
    bgPix = scribbled_sp_all_ind.bg ;
    fgPix = scribbled_sp_all_ind.fg ;
    ncenters = param.ncenter;
    GMMparam.ncenters = ncenters;
    X = A(:,bgPix);
    [Unaries, Mu, Sigma, ncentersB] = EM_init_kmeans(X, ncenters);
    [GMMparam.unarypotB, GMMparam.muB, GMMparam.sigmaB] = EM(X, Unaries, Mu, Sigma);
    GMMparam.ncentersB = ncentersB;

    X = A(:,fgPix);
    [Unaries, Mu, Sigma, ncentersF] = EM_init_kmeans(X, ncenters);
    [GMMparam.unarypotF, GMMparam.muF, GMMparam.sigmaF] = EM(X, Unaries, Mu, Sigma);
    GMMparam.ncentersF = ncentersF;
    labelParam = GMMparam;
    labelParam.fea_style = 'features';

    % initialization
    S = 0;
    start = zeros(n_img,1);
    stop = zeros(n_img,1);
    start(1)=1;
    
    
    d = size(histSP{1},1);
    h_bar = zeros(d,1);

    E0 = 0;
    e = cell(n_img,1);
    YS = cell(n_img,1);
    M_11 = cell(n_img,1);
    
    label_index = cell(n_img,1);
    Mglb = cell(n_img,1);
    Y_first = cell(n_img,1);
    
    
    % main part
    for i = 1:n_img
        feateVec_i = featureVector{i}';%d x n
        n_sp(i)=size(histSP{i},2);
        e{i} = ones(n_sp(i),1);
        
        % initialize YS with 0.6
        YS{i} = e{i}*0.6;

        M_11{i} = sparse(histSP{i}'*histSP{i});
        stop(i) = start(i) + n_sp(i) - 1;
        
        if i == n_img
            n = stop(i);
        else
            start(i+1) = stop(i) + 1;
        end  

        if i<=n_scri
            [ label_global_ind GMMprob] = get_labels_from_scribbles(feateVec_i,scribbled_sp_feat_all, scribbled_sp_all_ind,labelParam);
 
            label_local_ind = superpixel_labelInd{i};
            % keep global labels and local labels consistent
            logical_fg = zeros(n_sp(i),1);              logical_bg = zeros(n_sp(i),1);
            logical_fg(label_global_ind.fg) =1;         logical_bg(label_local_ind.bg) =1;
            temp_log = logical_fg+logical_bg;
            logical_fg(temp_log==2) = 0; 
            label_global_ind.fg = find(logical_fg);

            logical_fg = zeros(n_sp(i),1);              logical_bg = zeros(n_sp(i),1);
            logical_fg(label_local_ind.fg) = 1;         logical_bg(label_global_ind.bg) =1;
            temp_log = logical_fg + logical_bg;
            logical_bg(temp_log==2) = 0; 
            label_global_ind.bg = find(logical_bg);  

            %combine global labels and local labels 
            label_ind.fg = unique([label_local_ind.fg;label_global_ind.fg]);
            label_ind.bg = unique([label_local_ind.bg;label_global_ind.bg]);        

            [ Mglb{i} ,Y0{i}] = construct_intra_matrix( featureVector{i}', label_ind, edges_s{i} );

            Y1 = zeros(size(Y0{i})); Y1(label_global_ind.fg) = lambda0-lambda1;
            Mglb{i} = Mglb{i} - diag(Y1);  Y0{i} = Y0{i} - Y1;

             %GMM constraint
            label_index{i} = label_ind;

            % Add ECNN (Ejoint) energy
            maxVal=max(Mglb{i}\he{i}');   
            A = ((Mglb{i}\he{i}'))/(maxVal); 
            B = Mglb{i}\Y0{i};%+A;
            YS{i}=B;
            YS{i} = double(YS{i}>=0.5);
            h_bar = h_bar+histSP{i}*YS{i};

           %%%compute energy
            Mi = Mglb{i}-diag(Y0{i})+lambda2*M_11{i};
            ys_new = YS{i};
            Vi = lambda2*histSP{i}'*histSP{i}*ys_new;

            %%add ECNN
            ECNN=lambda3*he{i}*(1-ys_new);
            E0 = E0+ys_new'*Mi*ys_new-2*ys_new'*Vi+ECNN;

        else
            [ label_global_ind  GMMprob] = get_labels_from_scribbles( feateVec_i,scribbled_sp_feat_all, scribbled_sp_all_ind,labelParam );   
            [ Mglb{i} ,Y0{i}] = construct_intra_matrix( featureVector{i}', label_global_ind, edges_s{i} );
            Y1 = zeros(size(Y0{i})); Y1(label_global_ind.fg) = lambda0-lambda1;
            Mglb{i} = Mglb{i} - diag(Y1);  Y0{i} = Y0{i} - Y1;
             %GMM constraint
            label_index{i} = label_global_ind;  

             maxVal=max(Mglb{i}\he{i}');
            YS{i} = Mglb{i}\Y0{i}-(Mglb{i}\he{i}')/maxVal;          
            YS{i} = double(YS{i}>=0.5);     
            h_bar = h_bar+histSP{i}*YS{i};
            %%%compute energy
            Mi = Mglb{i}-diag(Y0{i})+lambda2*M_11{i};
            ys_new = YS{i};
            Vi = lambda2*histSP{i}'*histSP{i}*ys_new;

             %%add ECNN
            ECNN=lambda3*he{i}*(1-ys_new);
            E0 = E0+ys_new'*Mi*ys_new-2*ys_new'*Vi+ECNN;
 
        end
     end


    %% iterative part
    h_bar = h_bar/n_img;
    count = 1;
    options_c=optimset('Algorithm','interior-point','display','off');
    niter = 20;
    while count<=niter 

        E1 = 0;
        for i=1:n_img
            Mi = Mglb{i}-diag(Y0{i})+lambda2*M_11{i};
            Vi = lambda2*histSP{i}'*h_bar;
            ys_old = YS{i};        

            label_ind = label_index{i};
            lb = zeros(n_sp(i),1);lb(label_ind.fg)=1;
            ub = ones(n_sp(i),1);ub(label_ind.bg)=0;
            [ys_new,Ei] = fmincon(@(ye_old)myfun22(ys_old,Mi,Vi,he{i}),ys_old,[],[],[],[],lb,ub,[],options_c);

            Yp{i} = ys_new; 
            ys_new = double(ys_new>=0.5);   
            YS{i} =ys_new;
            h_bar = h_bar+histSP{i}*(ys_new-ys_old)/n_img;
            Vi = lambda2*histSP{i}'*h_bar;
            
            %%add ECNN
            ECNN=lambda3*he{i}*(1-ys_new);
            E1 = E1+ys_new'*Mi*ys_new-2*ys_new'*Vi+ECNN;
        end
        EE{count} = E1;

        fprintf('enger=%d\n',E1);
        count  = count+1;
        if abs(E1-E0)<0.001
            break
        else
            E0 = E1;
        end
    end
    fprintf('niter=%d\n',count-1);

    
    for i = 1:n_img
        YS{i} = double(YS{i}>=0.5);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%Calculate accuracy and jaccard similarity%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:n_img
        YS_temp = YS{i};
        [h,w,d] = size(labels{i});
        ind = find(YS_temp);
        Y_first = zeros(h,w);
        for j = 1:length(ind)
            Y_first(labels{i} == ind(j))=1;
        end
        Y_mask = uint8(repmat(Y_first,[1,1,3]));
        im = imread([img_path ,img_names{i},'.' imgstyle]);
        im = imresize(im,[224 224],'nearest');

       [imgMasks,segOutline,imgMarkup]=segoutput(im2double(im),double(Y_first+1));
    
        %% iCoseg
    %     groundtruth_path = [img_path,'GroundTruth/'];
    %     gtImage = imread([groundtruth_path,img_names{i},'.png']);
    %     groundtruth = double(gtImage(:,:,1)>0); groundtruth=imresize(groundtruth,[224 224],'nearest'); %unique(groundtruth); imshow(groundtruth);
    %     P(i) =sum(groundtruth(:)==Y(:)) ./ prod(size(groundtruth));
    %     Jar(i) =sum( (Y(:)==1) & (groundtruth(:)==1) ) ./ sum( (Y(:) | groundtruth(:))==1 );
    
        %% MSRC/Davis/SegTrackV2
        groundtruth_path = Gr_path;
        gtImage = imread([groundtruth_path,img_names{i},'.png']);
        gtImage(gtImage>200)=255; gtImage(gtImage<=200)=0;      
        groundtruth = double(gtImage(:,:,1))./255; unique(groundtruth);
        groundtruth=imresize(groundtruth,[224 224],'nearest'); %unique(groundtruth); imshow(groundtruth);
        P(i) =sum(groundtruth(:)==Y_first(:)) ./ prod(size(groundtruth));
          if sum((groundtruth(:)==1) )
            Jar(i) =sum( (Y_first(:)==1) & (groundtruth(:)==1) ) ./ sum( (Y_first(:) | groundtruth(:))==1 );
        else
            if sum((Y_first(:)==1))
                Jar(i) = 0;
            else
                Jar(i) =1;
            end            
        end
    %     Jar(i) =sum( (Y(:)==1) & (groundtruth(:)==1) ) ./ sum( (Y(:) | groundtruth(:))==1 );

        %% MFC
    %     groundtruth_path = [img_path,'GroundTruth\'];
    %     gtImage = imread([groundtruth_path,img_names{i},'.png']);
    %     groundtruth = double(gtImage>0 & gtImage<255);
    %     P(i) =sum(groundtruth(:)==Y(:)) ./ prod(size(groundtruth));
    %     Jar(i) =sum( (Y(:)==1) & (groundtruth(:)==1) ) ./ sum( (Y(:) | groundtruth(:))==1 );
    end
        MP = mean(P); MJ = mean(Jar);
        res(counter,:) = [MP,MJ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%NEW PART%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    superpixel_labelIndNew = find_best_alpha(superpixel_labelInd,YS,featureVector,param);
    
    AA = [];
    BB = [];
    for i=1:size(superpixel_labelIndNew,1)
        if ~isequal(superpixel_labelIndNew{i}.fg,superpixel_labelInd{i}.fg)
            AA = [AA,i];
        end
        if ~isequal(superpixel_labelIndNew{i}.bg,superpixel_labelInd{i}.bg)
            BB = [BB,i];
        end
    
    end

    superpixel_labelInd = superpixel_labelIndNew;
end