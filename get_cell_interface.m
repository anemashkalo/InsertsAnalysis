function [sorted_pluri] = get_cell_interface(pos_data,img_n,pxl_to_micron,chan_tmp,q_dist,norm)
% get the border between the cell types based on only one cell type being
% labelled
% still need dapi channel in all cells
clear mask;
i = img_n;% image number
mask = bwareafilt(pos_data(i).mask,[70 3000]);
% get the bw dist for the image
I_tmp = bwdist(mask);
%imshow(I_tmp,[]);
dat = pos_data(i).expression_data;
% then just bin the pixels from the processed masks to the bwdist output ?
pl_cell = pos_data(i).pluricells;
all_dist_pluricells=[];
indx = [];
all_dist_pluri = [];
norm_expr = [];
for jj=1:size(pl_cell,1)
indx(jj)=sub2ind(size(mask),round(dat(chan_tmp(1)).coord(pl_cell(jj),2)),round(dat(chan_tmp(1)).coord(pl_cell(jj),1)));% get the linear index of the pluri cells
% then find what distances these indexes correspond to
all_dist_pluricells(jj,1) = I_tmp((indx(jj)));% this is the distance that this indx(jj) pixel is away from the nearest prediff cell
if norm
norm_expr(jj,1) = dat(chan_tmp(q_dist)).int(pl_cell(jj))./dat(chan_tmp(1)).int(pl_cell(jj));
else
 norm_expr(jj,1) = dat(chan_tmp(q_dist)).int(pl_cell(jj));   
end
end
% assign distance to each pluri cell and combine data into matrix
sorted_pluri = [];
sorted_pluri = cat(2,dat(chan_tmp(1)).coord(pl_cell,1:2),all_dist_pluricells*pxl_to_micron,norm_expr);% last: normalized marker intensity

end



