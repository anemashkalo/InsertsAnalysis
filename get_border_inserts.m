function [boundary,slp,y_zero,fraction_wrong,prediff_are_above,outliers,border_esimate] = get_border_inserts(mask,mask_pluri,flag)
% this function works only if both cell labels are present, dose not use
% bwdist, but uses a linear fit through the border points
%close all
outliers=[];
mask_dil = imdilate(mask,strel('disk',2));
mask2_dil = imdilate(mask_pluri,strel('disk',2));
%imshowpair(mask_dil,mask2_dil);hold on
only_border =( mask_dil & mask2_dil);
%figure, imshow(only_border);
stats = regionprops(only_border,'Centroid','Area');
border_esimate = cat(1,stats.Centroid) ;
%hold on,plot(border_esimate(:,1),border_esimate(:,2),'.r')

%get rid of outliers
[r,~] = find(isoutlier(border_esimate,'gesd','ThresholdFactor',1)==1);%
outliers=border_esimate(r,:);
border_esimate(r,:)=[];

%hold on,plot(border_esimate(:,1),border_esimate(:,2),'*y')
% fit the line through these points
f=fitlm(border_esimate(:,1),border_esimate(:,2));
y_zero= f.Coefficients.Estimate(1);
slp = f.Coefficients.Estimate(2);

x_vect = 1:size(mask,1);
boundary = slp*x_vect+y_zero;
%hold on,plot(boundary);
%figure(1), hold on, plot(boundary,'b');

[fraction_wrong,prediff_are_above] = obtain_frac_inserts(mask,mask_pluri,slp,y_zero,flag);

end
