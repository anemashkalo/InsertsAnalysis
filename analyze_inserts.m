%% analyze the images from inserts

% TODO: generalize the code
% comment
% save variables
% automatic border detection

dir = 'D:\2019-04-16-Inserts_dapiBraTLC_Esi\DAPI-Bra-TLC-ESC\pos2';
chan = {'DAPI','BRA','TLC','ESI'};
%get raw images
raw_img=struct;
for jj=1:size(chan,2)
        raw_img(jj).dat =imread(['s_C00' num2str(jj) '.tif' ]);
end
% ff = dir('D:\2019-04-16-Inserts_dapiBraTLC_Esi\DAPI-Bra-TLC-ESC\pos5');
close all
k = 3;
prob_thresh=0.9;
ilastik_fn = ['s_C00' num2str(k) '_Probabilities.h5' ];
mask = readIlastikProbMask(ilastik_fn,prob_thresh);
k = 4;
ilastik_fn = ['s_C00' num2str(k) '_Probabilities.h5' ];
mask2 = readIlastikProbMask(ilastik_fn,prob_thresh);

imshowpair(mask,mask2);
imshowpair(raw_img(3).dat,raw_img(4).dat)
% get the border line between cells , coordinates of the line
%  TODO: determine the border automatically (from combination of two masks)
%y = m*x+b;
% pt1 = [60 486];%159 587 (pos5)
% pt2 = [526 887];%%847 748 (pos5)
% slope = (pt1(2)-pt2(2))/(pt1(1)-pt2(1));
% intrsept = pt1(2)-slope*pt1(1);
% x_vect = 1:size(raw_img(jj).dat,1);
% border = slope*x_vect+intrsept;
% hold on, plot(border,'c.');
%% background subrtact, get cell centroids, get only pluri cells centroids
close all
chan_tmp = [1 2 3 4]; % dapi bra prediffmarker plurimarker
ilastik_fn = ['s_C00' num2str(chan_tmp(1)) '_Probabilities.h5' ];
mask_allnuc = readIlastikProbMask(ilastik_fn,prob_thresh);
pxl_to_micron = 0.617;
% get background-subtracted image in dapi and bra channels
img_fin = struct;
for j=1:size(chan_tmp,2)
 img_fin(j).dat = simplebg([],mask_allnuc,raw_img(chan_tmp(j)).dat);
 %figure(j), imshow(img_fin(j).dat,[]);
end

expression_data = struct;
stats_tmp=[];
for i=1:size(chan_tmp,2)
stats_tmp = regionprops(mask_allnuc,img_fin(i).dat,'MeanIntensity','Centroid');
expression_data(i).coord = cat(1,stats_tmp.Centroid);
expression_data(i).int = cat(1,stats_tmp.MeanIntensity);

end
figure(2),imshow(raw_img(1).dat,[]); hold on%raw_img(2).dat,[]
%imshowpair(raw_img(1).dat,raw_img(2).dat);hold on

%  find the cells that are pluri cells
pluri_marker_thresh = 200; % 800 (pos 5) threshold for expression of YFP, to identify pluri cells
[pluricells,~]=find(expression_data(4).int > pluri_marker_thresh);
plot(expression_data(1).coord(pluricells,1),expression_data(1).coord(pluricells,2),'r.');
hold on
plot(border,'g');
size(pluricells)

%% get the cells at specific distances from the border
% strategy: move the line and see if any cell centroids intersect it
i = 1;
tmp = [];
cells_at_border = struct;
border = slope*x_vect+intrsept;
line_coord = cat(2,x_vect',border');
plot(line_coord(:,1),line_coord(:,2),'b');
for jj=1:size(pluricells,1)    
tmp= ipdm(expression_data(i).coord(pluricells(jj),1:2),line_coord,'Result','Structure','Subset','NearestNeighbor');% 
cells_at_border(jj).dist = tmp.distance;
cells_at_border(jj).bra = expression_data(2).int(pluricells(jj))/expression_data(1).int(pluricells(jj));% bra intensity normalized to dapi
plot(expression_data(i).coord(pluricells(jj),1),expression_data(i).coord(pluricells(jj),2),'m.');hold on

end

figure,plot(cat(1,cells_at_border.dist)*pxl_to_micron,cat(1,cells_at_border.bra),'b.');
ylabel('(BRA/DAPI) in pluri cells')
xlabel('distance from border,um')
title('Position 2')









