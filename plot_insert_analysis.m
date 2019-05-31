%% plot data from inserts
close all
matfile = 'D:\2019-05-22_inserts_noggin_IWP2\IWP2max1\IWP2_Inserts_BRAstain.mat';
load(matfile);
%% plot data from all images, bin the distances
%close all
bin_sz =10; % in microns
cell_dat = cat(1,pos_data.marker_vs_dist);% all pluri cells from all images
max_d = max(cell_dat(:,3));% in microns
bins = round(max_d/bin_sz);
bin_vect = (1:bin_sz:max_d);
bin_vect_fin = cat(2,bin_vect,max_d);
binned_data = [];
for jj=1:size(bin_vect_fin,2)-1
    tmp1 = [];
    tmp2 = [];
    %disp([bin_vect_fin(jj) bin_vect_fin(jj+1)])
    tmp1= find(cell_dat(:,3) <= bin_vect_fin(jj+1));
    tmp2 = find(cell_dat(:,3) > bin_vect_fin(jj));    
    binned_data(jj,1) = bin_vect_fin(jj);
    binned_data(jj,2) = mean(cell_dat(intersect(tmp1,tmp2),4));
    binned_data(jj,3) = std(cell_dat(intersect(tmp1,tmp2),4));
    binned_data(jj,4)= size(cell_dat(intersect(tmp1,tmp2),4),1);% how many cells contributed to the mean  
end
% plot data
 hold on,figure(1),errorbar(binned_data(:,1), binned_data(:,2), binned_data(:,3),'-*g','LineWidth',1.5);hold on
 %plot(binned_data(:,1), binned_data(:,2),'-*g','LineWidth',1.5);hold on
 h = figure(1);
     ylabel([ chan_nm{2} ' / DAPI  in pluri cells']);
     xlabel('Distance from cell types interface,um')
     title('All images');
     ylim([0 max(binned_data(:,2)+binned_data(:,3))]);     
h.CurrentAxes.LineWidth = 2;
h.CurrentAxes.FontSize = 14;

%% check the  cell distances that were determined
% load the mat file with the data
close all
q = 2; % image number (conting only the ones that were segmented in ilastik)
figure(q),imshowpair(pos_data(q).mask,pos_data(q).maskall);hold on
plot(pos_data(q).expression_data(chan_tmp(3)).coord(:,1),pos_data(q).expression_data(chan_tmp(3)).coord(:,2),'.r');hold on
plot(pos_data(q).marker_vs_dist(:,1),pos_data(q).marker_vs_dist(:,2),'c.');
text(pos_data(q).marker_vs_dist(:,1),pos_data(q).marker_vs_dist(:,2),num2str(pos_data(q).marker_vs_dist(:,3),2),'color','b');%all_dist_pluricell
title(num2str(pos_data(q).imgname));
