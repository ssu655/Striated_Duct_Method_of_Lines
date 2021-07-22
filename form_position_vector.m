function [CellPos, I, IntPos] = form_position_vector(seg_prop)
% This function extracts the centroids of cells and the positions of lumen
% segments and combines them in one position vector
%
% This function is solely for plotting purposes.

n_seg = size(seg_prop,2);
CellPos = cell(1, n_seg);
IntPos = cell(1, n_seg);

for j = 1:n_seg
    cell_prop = seg_prop{j}.cell_prop;
    lumen_prop = seg_prop{j}.lumen_prop;
    
    n_c = length(cell_prop);
    CellPos{j} = zeros(1,n_c);
    for i = 1:n_c
        % take the z coordinate of cell's centroid
        CellPos{j}(i) = cell_prop{i}.centroid(3);
    end
    
    IntPos{j} = lumen_prop.segment(1:end-1);

end
CellPos = cat(2,CellPos{:});
IntPos = cat(2,IntPos{:});

% since cells are not ordered along z axis, it needs to be sorted.
[CellPos,I] = sort(CellPos);