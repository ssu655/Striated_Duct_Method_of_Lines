function x = setup_IC(Conc, cell_prop, lumen_prop)

n_c = length(cell_prop);
n_l = lumen_prop.n_int;

x_c = zeros(9, n_c);
x_l = zeros(6, n_l);

CIC = cell2mat(struct2cell(Conc.CIC)); % convert the cellular IC to an array
LIC = cell2mat(struct2cell(Conc.LIC)); % convert the lumenal IC to an array

x_c(4:9,:) = repmat(CIC, 1, n_c);
x_c(1,:) = -26.1257;
x_c(2,:) = -52.2513;
x_c(3,:) = 1000;

x_l = repmat(LIC, 1, n_l);

x = [x_c(:); x_l(:)];
end