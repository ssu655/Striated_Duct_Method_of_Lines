function seg_prop = process_meshes(P_list, L_int, segment_meshes)

% This function takes an array of mesh files and process them one by one.
% 
% seg_prop is a cell array, each element is a structure
%
% seg_struct - cell_prop
%            - lumen_prop
%            

n_seg = length(segment_meshes.meshtypes);
seg_prop = cell(1,n_seg);

for i = 1:n_seg
    mesh_file = segment_meshes.filenames(i);
    [cell_geom, lumen_geom] = get_mesh_info(L_int, mesh_file);
    
    P = P_list{i};
    [cell_prop, lumen_prop] = process_mesh_info(cell_geom, lumen_geom, P);
    
    % if not the first duct segment, record the end position for the next
    if i ~= 1
        lumen_prop.segment = lumen_prop.segment + end_position;
        
        n_c = length(cell_prop);
        for j = 1:n_c
            % add the end position of prev duct segment to the cell
            % centriods of the next segment
            cell_prop{j}.centroid(3) = cell_prop{j}.centroid(3) + end_position;
            
        end
    end
    end_position = lumen_prop.segment(end);
    
    seg_struct = struct;
    seg_struct.cell_prop = cell_prop;
    seg_struct.lumen_prop = lumen_prop;
    seg_prop{i} = seg_struct;
end