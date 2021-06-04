addpath('C:\Users\lingm\Dropbox\PhD\method_of_lines_mesh\geometry_processing_package_repo\io')
addpath('C:\Users\lingm\Dropbox\PhD\method_of_lines_mesh\geometry_processing_package_repo\data')
addpath('C:\Users\lingm\Dropbox\PhD\method_of_lines_mesh\geometry_processing_package_repo\algebra')
addpath('C:\Users\lingm\Dropbox\PhD\method_of_lines_mesh\geometry_processing_package_repo\graphics')
addpath('C:\Users\lingm\Dropbox\PhD\method_of_lines_mesh\geometry_processing_package_repo\topology')
addpath('C:\Users\lingm\Dropbox\PhD\method_of_lines_mesh\geometry_processing_package_repo\parameterization')

addpath('C:\Users\lingm\Dropbox\PhD\method_of_lines_mesh\Duct_Cells_new')

% addpath('/Users/ssu655/Dropbox/PhD/method_of_lines_mesh/geometry_processing_package_repo/io')
% addpath('/Users/ssu655/Dropbox/PhD/method_of_lines_mesh/geometry_processing_package_repo/data')
% addpath('/Users/ssu655/Dropbox/PhD/method_of_lines_mesh/geometry_processing_package_repo/algebra')
% addpath('/Users/ssu655/Dropbox/PhD/method_of_lines_mesh/geometry_processing_package_repo/graphics')
% addpath('/Users/ssu655/Dropbox/PhD/method_of_lines_mesh/geometry_processing_package_repo/topology')
% addpath('/Users/ssu655/Dropbox/PhD/method_of_lines_mesh/geometry_processing_package_repo/parameterization')
% 
% addpath('/Users/ssu655/Dropbox/PhD/method_of_lines_mesh/Duct_Cells')

% ptCloud = pcread('PL-Cell_001.ply');
% pcshow(ptCloud)

filename = 'C:\Users\lingm\Dropbox\PhD\method_of_lines_mesh\mini_gland.ply';
[cells,faces,vertices] = read_ply_custom(filename);%(Files(1).name);
labels = unique(faces(:,4)); % array of unique face labels [api, bas, lat]
n_cell = size(cells,1);

figure 
hold on 

for i = 1:n_cell
    
    vert_start = cells(i,2)+1;
    vert_end = cells(i,2)+cells(i,1);
    face_start = cells(i,4)+1;
    face_end = cells(i,4)+cells(i,3);
    
    vertex = vertices([vert_start:vert_end],:);
    face = faces([face_start:face_end],:);
    face(:,1:3) = face(:,1:3) - vert_start + 1;
    
sf = size(face);   % face is a nf x 3 array
sv = size(vertex); % vertex is a nv x 3 array
labels = unique(face(:,4))

api_idx = find(face(:,4) == labels(1));
lat_idx = find(face(:,4) == labels(2));
bas_idx = find(face(:,4) == labels(3));
api_col = repmat([245,63,192]/255,[length(api_idx),1]);
lat_col = repmat([237,205,62]/255,[length(lat_idx),1]);
bas_col = repmat([71,230,111]/255,[length(bas_idx),1]);

col_arr = zeros(sf(1),3);
col_arr(api_idx,:) = api_col;
col_arr(bas_idx,:) = bas_col;
col_arr(lat_idx,:) = lat_col;

% write_obj('data/face.obj',face,vertex) % write mesh to obj format
po = patch('Faces',face(:,1:3),'Vertices',vertex,...
        'FaceColor','flat',...
        'FaceVertexCData',col_arr,...
        'CDataMapping','scaled');
view(45,20)
% 
[k, volume] = convhull(vertex);
% trisurf(k, vertex(:,1), vertex(:,2), vertex(:,3))

% identify the coordinate points of each face
face_coord = zeros(sf(1),9); %[x1, y1, z1, x2, y2, z2, x3, y3, z3]
face_area = zeros(sf(1),1);
for j = 1:sf(1)
    face_coord(j,:) = reshape(vertex(face(j,1:3),:)',1,[]);
    A = face_coord(j,1:3)-face_coord(j,4:6);
    B = face_coord(j,1:3)-face_coord(j,7:9);
    C = face_coord(j,4:6)-face_coord(j,7:9);
    face_area(j) = 0.5 * norm(cross(A,B));
end

api_area = sum(face_area(api_idx));
bas_area = sum(face_area(bas_idx));
lat_area = sum(face_area(lat_idx));

end

hold off 