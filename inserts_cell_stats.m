function [mask,mask_allnuc,img_fin,expression_data,pluricells,nofile] = inserts_cell_stats(datadir,img_pos,chan_tmp,prob_thresh,prediff_marker_thresh)
% gets all the stats from the suppled images and their ilastik probability
% masks; if image is not segmented with Ilasik, it is skipped
% images are background-subtracted here
% the pluri cells are identified as not prediff cells

 nofile = 0;
% mask for prediff cells
if ispc
ilastik_fn = [datadir '\' img_pos  '\s_C00' num2str(chan_tmp(3)) '_Probabilities.h5' ];%
else
   ilastik_fn = [datadir '/' img_pos  '/s_C00' num2str(chan_tmp(3)) '_Probabilities.h5' ]; %
end
disp(ilastik_fn);
if ~isfile(ilastik_fn)   
    disp('The image folder does not contain ilastik masks');
    mask1=[];
    mask_allnuc=[];
    img_fin=[];
    expression_data=[];
    pluricells=[];
    nofile = 1;
    return
    
else
mask1 = readIlastikProbMask(ilastik_fn,prob_thresh);
%mask = bwareafilt(mask1,[70 2000]);
% mask for all cells
if ispc
ilastik_fn = [datadir '\' img_pos '\s_C00' num2str(chan_tmp(1)) '_Probabilities.h5' ];
else
ilastik_fn = [datadir '/' img_pos '/s_C00' num2str(chan_tmp(1)) '_Probabilities.h5' ];    
end

mask_allnuc0 = readIlastikProbMask(ilastik_fn,prob_thresh);
mask_allnuc = bwareafilt(mask_allnuc0,[70 3000]);
%--------------- update mask prediff here only for this case (images
%may22,2019)% get rid of this section or generalize to cases where need to
%determine prediff mask from something else
 mask = [];
 mask2_tmp = [];
 mask2_tmp = imfill(mask1,'holes');
 mask3_tmp = [];
 mask3_tmp = imerode(mask2_tmp,strel('disk',5));
 mask0=[];
 mask0 = mask3_tmp&mask_allnuc;
 mask = bwareafilt(mask0,[70 3000]);
% imshowpair(mask3_tmp,mask );
%---------------------------
% get raw and background-subrtacted images for all channels
raw_img=struct;
img_fin = struct;
expression_data = struct;
stats_tmp=[];
for jj=1:size(chan_tmp,2)
        if ispc
        raw_img(jj).dat =imread([datadir '\' img_pos '\s_C00' num2str(jj) '.tif' ]);
        else
        raw_img(jj).dat =imread([datadir '/' img_pos '/s_C00' num2str(jj) '.tif' ]);    
        end
        
        img_fin(jj).dat = simplebg([],mask_allnuc,raw_img(chan_tmp(jj)).dat);
       % figure(jj), imshow(img_fin(jj).dat,[]);
        stats_tmp = regionprops(mask_allnuc,img_fin(jj).dat,'MeanIntensity','Centroid');% applied to background-subtracted images
        expression_data(jj).coord = cat(1,stats_tmp.Centroid);
        expression_data(jj).int = cat(1,stats_tmp.MeanIntensity);
end
%  find the cells that are pluri cells, as being non-diff cells
cell_id = [];
cell_id = (1:size(expression_data(1).coord,1))';% all cell ids;
cells_id_prediff = zeros(size(expression_data(1).coord,1),1);
[prediff_cells,~]=find(expression_data(chan_tmp(3)).int > prediff_marker_thresh);% find prediff cells (with marker above thresh)
cell_id(prediff_cells) = 0;
pluricells = nonzeros(cell_id);
% imshow(mask_allnuc,[]);
% hold on
% plot(expression_data(1).coord(:,1),expression_data(1).coord(:,2),'c.')
%  plot(expression_data(1).coord(prediff_cells,1),expression_data(1).coord(prediff_cells,2),'r.')
%  hold on
%  plot(expression_data(1).coord(pluricells ,1),expression_data(1).coord(pluricells ,2),'g.')
end
end