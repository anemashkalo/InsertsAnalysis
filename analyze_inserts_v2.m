%% get raw images, masks and border between cell types


close all
% input parameters:
chan_tmp = [1 2 3 4]; % [dapi stain prediff marker plurimarker], list channels in this order
prob_thresh=0.5; % make sure to segment dapi and prediff cells
prediff_marker_thresh = 150; %  threshold for expression of prediff nuc marker, to identify pluri cells, should be the same for all images
%datadir ='D:\2019-05-22_inserts_noggin_IWP2\IWP2max1';
datadir ='D:\2019-05-22_inserts_noggin_IWP2\nogginmax1';
%datadir = 'D:\2019-05-22_inserts_noggin_IWP2\controlmax1';

cd(datadir)
ff = dir(datadir);
chan_nm = {'DAPI','BRA','TLC','ESI'};
pxl_to_micron = 0.621;%0.617
% here can loop over images
pos_data = struct;
showimg =1; % show images and identified boundary
q = 1;
save_nm = 'noggin_Inserts_BRAstain';% provide full path
q_dist = 2; % channel for which to quantify the distance dependence (indexed into chan_tmp )
norm = 1;
for i=1:size(ff,1)   
if (strfind(ff(i).name,'pos')== 1)
disp(['processing image from  ' num2str(ff(i).name) ]); 
[mask,mask_allnuc,img_fin,expression_data,pluricells,nofile] = inserts_cell_stats(datadir,ff(i).name,chan_tmp,prob_thresh,prediff_marker_thresh);
if nofile == 1
    continue
end
pos_data(q).imgname = ff(i).name;
pos_data(q).maskall = mask_allnuc;
pos_data(q).mask = mask;% prediff cells
pos_data(q).expression_data =expression_data;
pos_data(q).pluricells = pluricells;% indexed into all cells
% get pluri cells' mean expression in channel (whatever was stained for) at specific distances from the border
[sorted_pluri] = get_cell_interface(pos_data,q,pxl_to_micron,chan_tmp,q_dist,norm);
[indx_tmp,~ ]= find(sorted_pluri(:,3)==0);% get rid of cells that were found to be zero microns away from interface
sorted_pluri(indx_tmp,:)=[];
%[cells_at_border] = marker_vs_border(chan_tmp,mask,slp,y_zero,expression_data,pluricells);
pos_data(q).marker_vs_dist = sorted_pluri;% distance already converted to microns
if showimg == 1
figure(q),imshowpair(pos_data(q).mask,pos_data(q).maskall);hold on
plot(pos_data(q).marker_vs_dist(:,1),pos_data(q).marker_vs_dist(:,2),'*');
text(pos_data(q).marker_vs_dist(:,1),pos_data(q).marker_vs_dist(:,2),num2str(pos_data(q).marker_vs_dist(:,3),2),'color','b');%all_dist_pluricell
title(num2str(pos_data(q).imgname));
figure(7+q),plot(pos_data(q).marker_vs_dist(:,3),pos_data(q).marker_vs_dist(:,4),'b.');
figure(7+q),title(num2str(pos_data(q).imgname));
xlabel('Distance from border, um');
xlabel([ 'Normalized ' (chan_nm(2)) 'expression' ]);
ylim([0 6]);
end
q = q+1;
end
end
% save the data for all images 

 %save(save_nm,'pos_data','chan_tmp','chan_nm','prob_thresh','prediff_marker_thresh','q_dist');
disp('done');

