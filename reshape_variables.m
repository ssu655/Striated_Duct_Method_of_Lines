function [x_c, x_l] = reshape_variables(x, seg_prop)
% This function reshape a long x vector to two matrix x_c and x_l
% 
% x is input variable of ODE system, shape: [1, N]
% N is total variable count
% 
% x_c [9, no. of cells in all duct segments]
% x_l [6, no. of lumen intervals in all ducts]


n_seg = size(seg_prop,2);
x_c = cell(1, n_seg);
x_l = cell(1, n_seg);

for j = 1:n_seg
    
    cell_prop = seg_prop{j}.cell_prop;
    lumen_prop = seg_prop{j}.lumen_prop;
    n_c = length(cell_prop);  % number of cells
    n_l = lumen_prop.n_int;   % number of lumen segments
    
    x_c{j} = reshape(x(1 : n_c*9),9,[]);        % [9, n_c]
    x_l{j} = reshape(x(1 + n_c*9 : n_c*9 + n_l*6),6,[]);  % [6, n_l]
    x = x(1 + n_c*9 + n_l*6:end);
    
end
x_c = cat(2,x_c{:});
x_l = cat(2,x_l{:});