function x = setup_IC(Conc, seg_prop)

% set up the initial condition to solve the system of ODE (f_ODE)

% x is [1 ; 9 * n_c + 6 * n_l] 
% 
% x_c(1,:) V_A
% x_c(2,:) V_B
% x_c(3,:) w_C
% x_c(4,:) Na_C
% x_c(5,:) K_C
% x_c(6,:) Cl_C
% x_c(7,:) HCO_C
% x_c(8,:) H_C
% x_c(9,:) CO_C
% 
% x_l(1,:) Na_A
% x_l(2,:) K_A
% x_l(3,:) Cl_A
% x_l(4,:) HCO_A
% x_l(5,:) H_A
% x_l(6,:) CO_A

n_seg = size(seg_prop,2);
x_seg = cell(1,n_seg);

for i = 1:n_seg
    cell_prop = seg_prop{i}.cell_prop;
    lumen_prop = seg_prop{i}.lumen_prop;
    n_c = length(cell_prop);  % number of cells
    n_l = lumen_prop.n_int;   % number of lumen segments

    x_c = zeros(9, n_c);
    x_l = zeros(6, n_l);

    CIC = cell2mat(struct2cell(Conc.CIC)); % convert the cellular IC to an array
    LIC = cell2mat(struct2cell(Conc.LIC)); % convert the lumenal IC to an array

    x_c(4:9,:) = repmat(CIC, 1, n_c);
    x_c(1,:) = -26.1257;
    x_c(2,:) = -52.2513;
    x_c(3,:) = 1000;

    x_l = repmat(LIC, 1, n_l);

    x_seg{i} = [x_c(:); x_l(:)];
end
x = cat(1,x_seg{:});
end