function [cell_geom, lumen_geom] = get_mesh_info(L)

% cell_geom is a cell array of struct containing the following features
% cell_raw    - api_idx      [1, n_api_mesh]
%             - bas_idx      [1, n_bas_mesh]
%             - lat_idx      [1, n_lat_mesh]
%             - face_coord   [9, n_mesh]
%             - face_area    [1, n_mesh]
%             - volume
%
% lumen_geom  - segment      [1, n_int+1]
%             - L
%             - start
%             - end
%             - length
%             - n_int
%             - radius
%             - volume


% add the path of the Duct_cells folder
addpath('/Users/ssu655/Dropbox/PhD/method_of_lines_mesh/Duct_Cells')
addpath('C:\Users\lingm\Dropbox\PhD\method_of_lines_mesh\Duct_Cells')

Files = dir('Duct_Cells/*.ply'); % a struct with ply file info
n_cell = length(Files)-1;
% index_i = [28,32];
cell_geom = cell(1,n_cell);

[face,~] = read_ply_custom(Files(1).name);
labels = unique(face(:,4)); % array of unique face labels [api, bas, lat]

start_z = 10000;
end_z = 0;

for i = 1:n_cell 
    
    [face,vertex] = read_ply_custom(Files(i).name);
    
    sf = size(face);
    
    api_idx = find(face(:,4) == labels(1));
    lat_idx = find(face(:,4) == labels(2));
    bas_idx = find(face(:,4) == labels(3));
    
    [~, volume] = convhull(vertex);
    
    face_coord = zeros(sf(1),9); %[x1, y1, z1, x2, y2, z2, x3, y3, z3]
    face_area = zeros(sf(1),1);
    for j = 1:sf(1)
        face_coord(j,:) = reshape(vertex(face(j,1:3),:)',1,[]);
        A = face_coord(j,1:3)-face_coord(j,4:6);
        B = face_coord(j,1:3)-face_coord(j,7:9);
        % C = face_coord(j,4:6)-face_coord(j,7:9);
        face_area(j) = 0.5 * norm(cross(A,B));
    end
    
    % record the cell centroid
    x = mean(mean(face_coord(:,[1,4,7])));
    y = mean(mean(face_coord(:,[2,5,8])));
    z = mean(mean(face_coord(:,[3,6,9])));
    
    cell_raw.centroid = [x,y,z];
    cell_raw.api_idx = api_idx;
    cell_raw.bas_idx = bas_idx;
    cell_raw.lat_idx = lat_idx;
    cell_raw.face_coord = face_coord;
    cell_raw.face_area = face_area;
    cell_raw.volume = volume;
    
    cell_geom{i} = cell_raw;
    
    
    % record the min and max z coordinate of this cell
    min_z = min(min(face_coord(:,[3,6,9])));
    max_z = max(max(face_coord(:,[3,6,9])));
    
    % update z limits of the lumen
    if min_z < start_z
        start_z = min_z;
    end
    if max_z > end_z
        end_z = max_z;
    end
    
end

lumen_geom.L = L;
lumen_geom.start = start_z;
lumen_geom.end = end_z;
lumen_geom.length = end_z - start_z;
n_int = ceil(lumen_geom.length / L);
lumen_geom.n_int = n_int;

lumen_geom.segment = L*[0:(n_int)] + start_z;

x = cell_geom{end}.face_coord(cell_geom{end}.api_idx,[1,2]);
y = sqrt(x(:,1).^2 + x(:,2).^2);
z = floor(mean(y));
lumen_geom.radius = 4.0824829;%z; % um
lumen_geom.volume = pi*lumen_geom.radius^2*L; % um^3

end