function [cell_prop, lumen_prop] = process_mesh_info(cell_geom, lumen_geom, P)

% cell_prop is a cell array of struct containing the following features
% cell_struct - volume
%             - centriod      [1, 3]
%             - loc_int       [1, n_loc_int]
%             - api_area_int  [1, n_loc_int]
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
% lumen_prop - segment [1, n_int+1]
%            - volume
%            - X_area
%            - n_int
%            - L


n_cell = length(cell_geom);

cell_prop = cell(1,n_cell);
lumen_segment = lumen_geom.segment;

% the apical area used to scale the conductances G
A = 104.719755;


for i = 1:n_cell
    cell_raw = cell_geom{i};
    
    % calculating areas
    api_face_area = cell_raw.face_area(cell_raw.api_idx);
    api_area = sum(api_face_area);
    bas_area = sum(cell_raw.face_area(cell_raw.bas_idx));
    lat_area = sum(cell_raw.face_area(cell_raw.lat_idx));
    
    % calculate the mean z coordinate of each apical triangle
    api_coord_z = cell_raw.face_coord(cell_raw.api_idx,[3,6,9]);
    api_coord_z_mean = mean(api_coord_z,2);
    
    % sort apical triangles into corresponding lumen segments
    api_lumen_conn = discretize(api_coord_z_mean,lumen_segment);
    
    % find the unique lumen segments for this cell
    loc_int = unique(api_lumen_conn');
    api_area_int = zeros(size(loc_int));
    
    % loop through the lumen segments to calculate areas per segment
    for j = 1:length(loc_int)
        api_int_ind = find(api_lumen_conn == loc_int(j));
        api_area_int(j) = sum(api_face_area(api_int_ind));
    end
    %loc_int = 1:5;
    %api_area_int = 0.2*api_area*ones(size(loc_int));
    
    cell_struct.volume = cell_raw.volume;
    cell_struct.centroid = cell_raw.centroid;
    cell_struct.loc_int = loc_int;
    cell_struct.api_area_int = api_area_int;
    cell_struct.api_area = api_area;
    cell_struct.baslat_area = bas_area + lat_area;
    cell_struct.api_lumen_conn = api_lumen_conn;
    
    A_A = api_area;
    A_B = bas_area + lat_area;
    
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
    
    cell_prop{i} = cell_struct;
end

lumen_prop.segment = lumen_geom.segment;
lumen_prop.volume = pi * lumen_geom.radius^2 * lumen_geom.L; % um^3
lumen_prop.X_area = pi * lumen_geom.radius^2; % um^2
lumen_prop.n_int = lumen_geom.n_int;
lumen_prop.L = lumen_geom.L;

end