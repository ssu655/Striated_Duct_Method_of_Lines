function M = mass(t, x, cell_prop, lumen_prop)

% This function sets up the mass matrix M for the system of ODE
% M_c [9,9] is the mass matrix of one cell
% M_l [6,6] is the mass matrix of one lumenal segment
% M is arranged as follows:
% M = | M_c                            |
%     |     M_c                        |
%     |         M_c                    |
%     |            ...                 |
%     |                M_l             |
%     |                    M_l         |
%     |                        M_l     |
%     |                            ... |
%       _____________  ______________
%            n_c              n_l


n_c = length(cell_prop);
n_l = lumen_prop.n_int;

x_c = reshape(x(1 : n_c*9),9,[]); %[9, n_c]
% x_l = reshape(x(1+n_c*9 : end),6,[]); %[6, n_l]
% disp(x_c)
% disp(x_l)
% n_l = size(x_l,2);
    
m_l = eye(6);
M_l = repmat({m_l}, 1, n_l);

m_c = eye(9);
m_c(1,1) = 0;
m_c(2,2) = 0;
M_c = repmat({m_c}, 1, n_c);
V = 1:10:9*9; % pick the diagonal elements of m_c matrix
e = [4,5,6,7,8,9];

for i = 1:n_c
    M_c{i}(V(e)) = x_c(3,i); % place volume on diagonal
    M_c{i}(e,3) = x_c(e,i); % place ion concentration on corresponding places
end

M = [M_c,M_l];
M = blkdiag(M{:});

end