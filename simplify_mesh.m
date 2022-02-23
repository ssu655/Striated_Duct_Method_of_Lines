function [cell_prop, lumen_prop] = simplify_mesh(cell_prop, lumen_prop)

sim_lumen_prop = struct;

sim_lumen_prop.n_disc = 2;

% the lumen disc where we split SD and ID, SD is from node 0 to split node
lumen_split = 90;

% summing the SD length and ID length
disc_length = zeros(1,2); % [SD, ID]
disc_length(1) = sum(lumen_prop.disc_length(1:lumen_split));
disc_length(2) = sum(lumen_prop.disc_length(lumen_split+1:end));
sim_lumen_prop.disc_length = disc_length;

disc_volume = zeros(1,2); % [SD, ID]
disc_volume(1) = sum(lumen_prop.disc_volume(1:lumen_split));
disc_volume(2) = sum(lumen_prop.disc_volume(lumen_split+1:end));
sim_lumen_prop.disc_volume = disc_volume;

sim_lumen_prop.disc_X_area = disc_volume./disc_length;

% there are 3 nodes: 0,1,2; [0,1] is SD; [1,2] is ID
sim_lumen_prop.disc_out_Vec = [0, 1]; 

sim_cell_prop = cell(1, 2);

n_c = length(cell_prop); % [ID cell, SD cell]

% initiating the simplified cell_prop cell struct
for i = [1, 2]
    sim_cell_prop{i} = struct;
    sim_cell_prop{i}.api_area = 0;
    sim_cell_prop{i}.baslat_area = 0;
    sim_cell_prop{i}.api_area_discs = [0, 0];
end
sim_cell_prop{1}.type = 'I';
sim_cell_prop{2}.type = 'S';

n_S = 0;
n_I = 0;

for i = 1:n_c
    if cell_prop{i}.type == 'S'
        j = 2;
        n_S = n_S + 1;
    else
        j = 1;
        n_I = n_I + 1;
    end
    
    sim_cell_prop{j}.api_area = sim_cell_prop{j}.api_area + cell_prop{i}.api_area;
    sim_cell_prop{j}.baslat_area = sim_cell_prop{j}.baslat_area + cell_prop{i}.baslat_area;
    sim_cell_prop{j}.api_area_discs(1) = sim_cell_prop{j}.api_area_discs(1) + sum(cell_prop{i}.api_area_discs(1:lumen_split));
    sim_cell_prop{j}.api_area_discs(2) = sim_cell_prop{j}.api_area_discs(2) + sum(cell_prop{i}.api_area_discs(lumen_split+1:end));
end

% mean distance of cell from duct node 0
sim_cell_prop{1}.mean_dist = 100; % ID cell
sim_cell_prop{2}.mean_dist = 50; % SD cell

% record how many cells comprise the sim_cell
sim_cell_prop{1}.n_c = n_I; % ID cell
sim_cell_prop{2}.n_c = n_S; % SD cell

cell_prop = sim_cell_prop;
lumen_prop = sim_lumen_prop;
end