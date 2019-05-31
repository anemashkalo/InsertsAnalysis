function [cells_at_border] = marker_vs_border(chan_tmp,mask,slp,y_zero,expression_data,pluricells)
% used only with the data, where the cell interface is represented by a
% line, identified based on marker from both cell types (pluri and prediff
% cells are both labelled)
x_vect =(1:size(mask,1));
tmp = [];
cells_at_border = struct;
interface = slp*x_vect+y_zero;
line_coord = cat(2,x_vect',interface');
%plot(line_coord(:,1),line_coord(:,2),'b');
for jj=1:size(pluricells,1)    
tmp= ipdm(expression_data(chan_tmp(1)).coord(pluricells(jj),1:2),line_coord,'Result','Structure','Subset','NearestNeighbor');% 
cells_at_border(jj).dist = tmp.distance;
cells_at_border(jj).bra = expression_data(chan_tmp(2)).int(pluricells(jj))/expression_data(chan_tmp(1)).int(pluricells(jj));% bra intensity normalized to dapi
%plot(expression_data(chan_tmp(1)).coord(pluricells(jj),1),expression_data(chan_tmp(1)).coord(pluricells(jj),2),'m.');hold on

end





end