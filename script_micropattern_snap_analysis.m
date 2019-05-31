%% analyse sorted pattern separate time points using the inserts analysis code
% for time points, where the cell border is established, get the Smad4 as a
% function of distance to cell border, 

% get all time point masks for each channel 
tic
simple_seg = 0;
 dir_img ='E:\allSortingData\2017-09-25-liveSortingwithRegistration\MIPs_1st20hrs';% E:\allSortingData\2017-07-14-Smad4sorting_maxProjections
 pos =1;
 tgroup = [];
 chan1=1;
 foreground = 1;%1
 thresh=0.5; 
 [mask_prediff,reader1] = get_mask(dir_img,chan1,pos,tgroup,thresh,foreground,simple_seg); 
 chan1=2; 
 foreground = 1;%1
 [mask_pluri,reader2] = get_mask(dir_img,chan1,pos,tgroup,thresh,foreground,simple_seg);
 foreground = 1;
  chan1=3; 
 [~,reader3] = get_mask(dir_img,chan1,pos,tgroup,thresh,foreground,simple_seg);
  toc 
     nT = size(mask_prediff,3)  
  
     %% get raw images and apply bwdist
  % TODO: get nuc to cyto, clean the plui cells by size
    % close all
  pos_data= struct;
  nz = 1;
  tp = 80;  
  dil = 5;
  mask_pluri(:,:,tp)=bwareafilt(mask_pluri(:,:,tp),[100 500]);
  figure,imshowpair(mask_prediff(:,:,tp),mask_pluri(:,:,tp));
  chan1 = 1;
  iPlane1=reader1.getIndex(nz - 1,chan1 -1 , tp - 1) + 1;%chan1 -1
  img_raw_prediff=bfGetPlane(reader1,iPlane1);  
  iPlane2=reader2.getIndex(nz - 1,chan1 -1 , tp - 1) + 1;%chan1 -1
  img_raw_pluri=bfGetPlane(reader2,iPlane2);
 % figure, imshowpair(img_raw_prediff,img_raw_pluri);  
  iPlane3=reader3.getIndex(nz - 1,chan1 -1 , tp - 1) + 1;%chan1 -1
  img_raw_var=bfGetPlane(reader3,iPlane3);
  %figure, imshowpair(img_raw_prediff,img_raw_var);  
  
   mask_cyto1 = imdilate(mask_pluri(:,:,tp),strel('disk',dil));
   mask_cyto = mask_cyto1 &~mask_pluri(:,:,tp);
%   figure, imshow(mask_cyto,[]);
  
%%
%close all
  q = 1;
% todo: clean up the junk; get the nuc/cyto SMAD4
pos_data(q).imgname = tp;
pos_data(q).maskall = (mask_prediff(:,:,tp) + mask_pluri(:,:,tp));
pos_data(q).mask = mask_prediff(:,:,tp);% prediff cells
pos_data(q).mask_pluri = mask_pluri(:,:,tp);% 

%chan_tmp = [plurimarker smad4 prediffmarker];
raw_img = struct;
raw_img(1).dat = img_raw_pluri;
raw_img(2).dat = img_raw_var;
raw_img(3).dat = img_raw_prediff;

img_fin=struct;
expression_data=struct;
 % [dapi(norm.marker) stain prediff marker plurimarker], list channels in this order
chan_tmp = [1 2 3];%2 3 1
for jj=1:size(chan_tmp,2)
        img_fin(jj).dat = simplebg([],pos_data(q).mask_pluri,raw_img(jj).dat);        
        figure(jj), imshow(raw_img(jj).dat,[]);
        stats_tmp = regionprops(pos_data(q).mask_pluri,raw_img(jj).dat,'MeanIntensity','Centroid','PixelIdxList');% applied to background-subtracted images
        expression_data(jj).coord = cat(1,stats_tmp.Centroid);
        expression_data(jj).int = cat(1,stats_tmp.MeanIntensity);
        if jj == 2
            for k=1:size(stats_tmp,1)
                tmp_mask = zeros(size(mask_pluri(:,:,tp)));
                tmp_mask(stats_tmp(k).PixelIdxList) = 1;
                tmp_mask2 = imdilate(tmp_mask,strel('disk',dil));                
                mask_cyto3 = tmp_mask2 &~tmp_mask;
                %figure,imshowpair(mask_cyto3,tmp_mask);
                
        stats_tmp1 = regionprops(mask_cyto3,raw_img(jj).dat,'MeanIntensity');%
        expression_data(jj).S4cyto(k,1) = stats_tmp1.MeanIntensity;
        stats_tmp2 = regionprops(tmp_mask,raw_img(jj).dat,'MeanIntensity');%
        expression_data(jj).S4nuc(k,1) = stats_tmp2.MeanIntensity;
            end
        end
end
chan_tmp = [1 2 3 4];%2 3 1
norm = 0;% normalize or not the channel for which to quantify distance
s4chan_bg = 1100;%1100
expression_data(4).int = (expression_data(2).S4nuc -s4chan_bg)./(expression_data(2).S4cyto-s4chan_bg);
pos_data(q).expression_data =expression_data;
pos_data(q).pluricells = (1:size(expression_data(1).coord,1))';% indexed into all cells
% get pluri cells' mean expression in channel (whatever was stained for) at specific distances from the border
chan_nm = {'RFPplurimarker','nucSMAD4','CFPprediff'};
pxl_to_micron = 0.617;

%save_nm = 'Inserts_BRAstain';% provide full path
q_dist =4; % channel for which to quantify the distance dependence (indexed into chan_tmp )
[sorted_pluri] = get_cell_interface(pos_data,q,pxl_to_micron,chan_tmp,q_dist,norm);
%[cells_at_border] = marker_vs_border(chan_tmp,mask,slp,y_zero,expression_data,pluricells);

pos_data(q).marker_vs_dist = sorted_pluri;% distance already converted to microns

figure(q),imshowpair(mask_cyto,pos_data(q).mask_pluri);hold on%pos_data(q).mask
plot(pos_data(q).marker_vs_dist(:,1),pos_data(q).marker_vs_dist(:,2),'*');
text(pos_data(q).marker_vs_dist(:,1),pos_data(q).marker_vs_dist(:,2),num2str(pos_data(q).marker_vs_dist(:,3),2),'color','b');%all_dist_pluricell
title(num2str(pos_data(q).imgname));
figure(7+q),plot(pos_data(q).marker_vs_dist(:,3),pos_data(q).marker_vs_dist(:,4),'.');hold on
figure(7+q),title(num2str(pos_data(q).imgname));
xlabel('Distance from border, um');
%ylabel([ 'Normalized ' (chan_nm(2)) 'expression' ]);
ylim([0.5 1.5]);




