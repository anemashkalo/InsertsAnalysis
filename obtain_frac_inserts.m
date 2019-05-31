function [fraction_wrong,prediff_are_above] = obtain_frac_inserts(mask,mask_pluri,slp,y_zero,flag)


% determine which way are the prediff cells (below or above the line) and
% fraction wrong

fraction_wrong = [];
prediff_are_above = [];

q = 1;
prediff_stats = regionprops(mask,'Centroid');
all_cells_prediff = cat(1,prediff_stats.Centroid);
cells_above_border=[];
cells_below_border=[];
test_pt=[];
for jj=1:size(all_cells_prediff,1)
    test_pt = slp*all_cells_prediff(jj,1)+y_zero;
    if all_cells_prediff(jj,2) < test_pt
        cells_above_border(q,1:2) = all_cells_prediff(jj,1:2);
        q = q+1;
    end
end



fraction_wrong = size(cells_above_border,1)/size(all_cells_prediff,1);
prediff_are_above = 0;
if fraction_wrong > 0.5
    prediff_are_above = 1;
end
if flag == 1
figure,imshowpair(mask,mask_pluri);hold on
plot(cells_above_border(:,1),cells_above_border(:,2),'bp');
plot(all_cells_prediff(:,1),all_cells_prediff(:,2),'r.');
end
end
