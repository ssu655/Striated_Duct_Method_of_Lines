function segment_meshes = gather_mesh_files()

segment_meshes = struct;

% form an array of mesh file names, one file for each duct segment
directory = 'C:\Users\lingm\Dropbox\PhD\method_of_lines_mesh\';
% directory = '/Users/ssu655/Dropbox/PhD/method_of_lines_mesh/';
segment_meshes.filenames = [ ...
    append(directory, "mini_gland_intercalated_1.ply") , ...
     append(directory, "mini_gland_slurm-20373309_1.ply"), ...
     append(directory, "mini_gland_slurm-20373309_1.ply") ...
     ];

% form an array of types corresponding to the duct segments.
% 1 is intercalated duct
% 2 is striated duct
segment_meshes.meshtypes = [1,2,2];

% form a matrix describing which duct segment feeds into which other segment
% M(i,j) = 1 means duct segment j feeds into segment i
% M(i,j) = 0 means segment j does not feed into segment i
% M(i,:) = 0 means segment i is connected to an acinus, ie primary saliva
% M matrix has to be a lower triangular matrix, i.e. 1 can only be found
% below the diagonal elements.
% segment_meshes.connectivity = [0,0;1,0];
end