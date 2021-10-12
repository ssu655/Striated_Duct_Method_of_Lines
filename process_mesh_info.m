function [S_cell_prop, lumen_prop] = process_mesh_info(L,P)

% S_cell_prop is a cell array of struct containing the following features
% cell_struct - mean_dist    (mean distance of cell from node 0)
%             - api_area_discs  [1, n_disc]
%             - api_area
%             - baslat_area
%             - scaled_rates  - L_A
%                             - L_B
%                             - G_ENaC
%                             - G_CFTR
%                             - G_BK
%                             - G_K_B
%                             - G_P_Na
%                             - G_P_K
%                             - G_P_Cl
%                             - NKA   - alpha_A 
%                                     - alpha_B
% lumen_prop - n_disc
%            - disc_length [1, n_disc]  (discs are index from node 0/duct root)
%            - disc_volume [1, n_disc]
%            - disc_X_area [1, n_disc]
%            - disc_out_Vec [1, n_disc]
%            - 
%%
D_fname = dir("*/*_duct.ply");
D_fname = strcat(D_fname.folder, filesep, D_fname.name);
[nodes, segments, radii, ~] = read_duct_tree_mesh(D_fname);

n_seg = size(segments,1); % segment is the sections of duct provided in the mesh file

% segments are indexed from node 0 or duct outlet
seg_length = vecnorm(nodes(segments(:,1)+1,:) - nodes(segments(:,2)+1,:),2,2);

% discs are further discretisation of the duct segment
n_disc = 0; 

% discs are also indexed from node 0
disc_length = zeros(1,6); % [1, n_disc] starting the vector a bit short, will get longer

% which segment the disc belongs to
d_s_Vec = zeros(1,6); % [1, n_disc]
disc_X_area = zeros(1,6); 

% keep track of the output segment/disc of each segment/disc, in terms of water flow
seg_out_Vec = zeros(1, n_seg);
disc_out_Vec = zeros(1, 6);

for i = 1:n_seg
    n = ceil(seg_length(i)/L); % num of discs in this segment
    
    % the first n discs has length L
    disc_length(n_disc+1 : n_disc+n-1 ) = L;
    
    % the last disc has the remainder as length
    if mod(seg_length(i), L)~= 0 
        disc_length(n_disc+n) = mod(seg_length(i), L);
    else
        disc_length(n_disc+n) = L;
    end
    
    % record the duct segment the discs belongs to
    d_s_Vec(n_disc+1 : n_disc+n) = i;
    
    % dis_mid_point is used to interpolate disc radius
    disc_mid_point = zeros(1,n);
    for j = 1:n
        disc_mid_point(j) = disc_length(j)/2+sum(disc_length(1:j-1));
    end
    disc_X_area(n_disc+1 : n_disc+n) = pi*(radii(i) + (radii(i+1)-radii(i))./seg_length(i).*disc_mid_point).^2;
    
    % seg_out [1, n_seg] 
    % seg_out(i) is the output segment of duct segment i
    % i.e. seg_out(2) = 1, the output segment of segment 2 is 1.
    seg_out = find(segments(:,1) == segments(i,2)); 
    % find the input node that equals the output node of the current segment.
    
    % if there is an output segment, ot is it not the last segment
    if seg_out ~= 0 
        seg_out_Vec(i) = seg_out;
        
        % find the disc index of the corresponding output segment 
        seg_out_disc = find(d_s_Vec == seg_out);
        % the output disc of the first dist is the last disc of the output segment
        % (first - close to node 0, last - far from node 0)
        disc_out_Vec(n_disc+1) = seg_out_disc(end);
    end
    % the other discs in the segment are consequtive.
    disc_out_Vec(n_disc+2:n_disc+n) = n_disc+1:(n_disc+n-1); 
    
    % update total num of discs with those in the current segment
    n_disc = n_disc + n;
    
end
disc_volume = disc_X_area .* disc_length;

% Saving function output
%   - duct segments information is not saved pass this function
%   - only duct discs information in returned as function output
lumen_prop.disc_length = disc_length;
lumen_prop.n_disc = n_disc;
lumen_prop.disc_volume = disc_volume;
lumen_prop.disc_X_area = disc_X_area;
lumen_prop.disc_out_Vec = disc_out_Vec;

%%
S_fnames = dir("*/*Cell_S*.ply");
n_S_cell = length(S_fnames);
S_cell_prop = cell(1,n_S_cell);

for i = 1:n_S_cell 
    S_fname = strcat(S_fnames(i).folder, filesep, S_fnames(i).name);
    
    cell_struct = process_cell_info(S_fname, seg_out_Vec, seg_length, d_s_Vec, disc_length);
     
    % the apical area used to scale the conductances G
    A = 104.719755;
    
    A_A = cell_struct.api_area;
    A_B = cell_struct.baslat_area;
    
    % scale the rates based on cell surface areas
    scaled_rates = struct;
    scaled_rates.L_A    = P.L_A * A / A_A;
    scaled_rates.L_B    = P.L_B * A / A_B;
    scaled_rates.G_ENaC = P.G_ENaC * A / A_A;
    scaled_rates.G_CFTR = P.G_CFTR * A / A_A;
    scaled_rates.G_BK   = P.G_BK * A / A_A;
    scaled_rates.G_K_B  = P.G_K_B * A / A_B;
    scaled_rates.G_P_Na = P.G_P_Na * A / A_A;
    scaled_rates.G_P_K  = P.G_P_K * A / A_A;
    scaled_rates.G_P_Cl = P.G_P_Cl * A / A_A;
    scaled_rates.NKA    = P.NKA;
    scaled_rates.NKA.alpha_A = P.NKA.alpha_A * A / A_A;
    scaled_rates.NKA.alpha_B = P.NKA.alpha_B * A / A_B;
    cell_struct.scaled_rates = scaled_rates;
    
    S_cell_prop{i} = cell_struct;
end

end

function cell_struct = process_cell_info(S_fname, seg_out_Vec, seg_length, d_s_Vec, disc_length)

[verts, faces, tri_type, duct_idx, ~, dist_along_duct] = read_striated_intercalated_mesh(S_fname);

% duct_idx - segment index for the triangle face starting from 0
duct_idx = duct_idx+1; % to match the zero based indexing of mesh file

% calculate individual triangle face area
n_faces = length(tri_type);
face_area = zeros(1,n_faces);

for j = 1:n_faces
    A = verts(faces(j,1),:) - verts(faces(j,2),:);
    B = verts(faces(j,2),:) - verts(faces(j,3),:);
    face_area(j) = 0.5 * norm(cross(A,B));
end

api_idx = find(tri_type == 0); % [1, n_apical_idx]
lat_idx = find(tri_type == 1); % [1, n_lateral_idx]
bas_idx = find(tri_type == 2); % [1, n_basal_idx]

api_binary = zeros(size(duct_idx)); % [1, n_triangle_idx]
api_binary(api_idx) = 1;
duct_idx_api = duct_idx.*api_binary; % [1, n_triangle_idx]

api_face_area = face_area.*api_binary'; % [1, n_triangle_idx]
api_area = sum(api_face_area);
bas_area = sum(face_area(bas_idx));
lat_area = sum(face_area(lat_idx));

% the duct index corresponding to the api triangles
duct_seg = unique(duct_idx(api_idx));

% Calculating the total apical triangle area corresponding to each disc
n_disc = length(disc_length);
api_area_discs = zeros(1,n_disc); 

% first find out the segments covered by the cell
duct_seg_tot = unique(duct_idx);

% convert distance along duct segment to distance from node 0
total_dist_along_duct = zeros(size(duct_idx));

% Since triangle dist_along_duct is segment based, we need to consider segment by segment
for j = 1:length(duct_seg_tot)
    s = duct_seg_tot(j);
    tri_seg_idx = find(duct_idx == s);
    dist_start_seg = calc_dist_start_seg(s, seg_out_Vec, seg_length);
    
    % calculate the total distance
    %       - discs are indexed from node 0 of the whole duct
    %       - in contrast, the triangle's distance along duct is measured from acinus
    %       - thus dist along duct needs to be reversed.
    total_dist_along_duct(tri_seg_idx) = seg_length(s) - dist_along_duct(tri_seg_idx) + dist_start_seg;
end

% the average distance along the duct for a cell
mean_dist = mean(total_dist_along_duct);

% For loop to put apical triangles into the corresponding discs, based on
% their total distance along duct.
% This needs to be done segment by segment, in case of duct branches.
% (triangles on two branches have similar distance but belong to diff discs)
for j = 1:length(duct_seg)
    
    s = duct_seg(j); 
    % apical indices corresponding to segment s
    tri_seg_idx = find(duct_idx_api == s);
    total_dist_along_duct_seg = total_dist_along_duct(tri_seg_idx);
    
    % find all discs in segment s
    all_discs_in_seg = find(d_s_Vec == s);
    distal_disc = min(all_discs_in_seg); % close to node 0 disc
    proxim_disc = max(all_discs_in_seg); % far from node 0 disc
    
    % make an array of bins for apical triangle distances, using the discs
    % consider all the discs from node 0 to segment s
    disc_edges = zeros(1,(proxim_disc + 1));
    dist_start_seg = calc_dist_start_seg(s, seg_out_Vec, seg_length);
    for k = distal_disc : proxim_disc + 1
        disc_edges(k) = sum(disc_length(distal_disc:k-1)) + dist_start_seg;
    end
    
    % drop the distance along duct into the dise edge bins
    api_disc_conn = discretize(total_dist_along_duct_seg, disc_edges);
    
    % array of the discs that has apical triangles in them
    loc_disc = unique(api_disc_conn);
    
    % the apical triangle areas corresponding the segment j
    api_face_area_seg = api_face_area(tri_seg_idx);
    
    % summing up all the triangle area in each bin/disc
    for k = 1:length(loc_disc)
        api_disc_ind = find(api_disc_conn == loc_disc(k));
        api_area_discs(loc_disc(k)) = sum(api_face_area_seg(api_disc_ind));
    end
end

cell_struct.mean_dist = mean_dist;
cell_struct.api_area_discs = api_area_discs;
cell_struct.api_area = api_area;
cell_struct.baslat_area = bas_area + lat_area;
end

function dist_start_seg = calc_dist_start_seg(seg_idx, seg_out_Vec, seg_length)
% distance_start = distance from node 0 at the start of seg_idx
dist_start_seg = 0;

% starting from the current segment, tracing its output segment until the root
k = seg_idx;
while seg_out_Vec(k) ~= 0
    k = seg_out_Vec(k);
    
    % adding up the output segment lengths
    dist_start_seg = dist_start_seg + seg_length(k);
end
end