%% read the dynamic simulation data z and post process it according to 
% visualisation requirement

% Cell = [Va Vb w_C Na_C K_C Cl_C HCO_C H_C CO_C]
% Lumen = [Na_L K_L Cl_L HCO_L H_L CO_L]
% z =[time_span; Cell * n_C + Lumen * n_l]
zz = z;
n_t = size(z,1);

% convert H_C, H_L to pH
n_c = length(s_cell_prop);
for i = 0:n_c-1
    H_C = z(:,8+9*i);
    pH_C = -log10(H_C*1e-3);
    zz(:,8+9*i) = pH_C;
end
n_l = s_lumen_prop.n_disc;
for i = 0:n_l-1
    H_L = z(:,n_c*9 + 5+6*i);
    pH_L = -log10(H_L*1e-3);
    zz(:,n_c*9 + 5+6*i) = pH_L;
end

% take out the Va, Vb, w_C, CO_C, CO_L information
zzz = zeros(n_t,10);
cell_remain = [4:8];
n_c_remain = length(cell_remain);
for i = 0:n_c-1
    remains = zz(:,cell_remain+9*i);
    zzz(:,[1:n_c_remain]+n_c_remain*i) = remains;
end

lumen_remain = [1:5];
n_l_remain = length(lumen_remain);
for i = 0:n_l-1
    remains = zz(:,n_c*9 + lumen_remain +6*i);
    zzz(:,[1:n_c_remain] + n_c*n_c_remain + n_l_remain*i) = remains;
end

% insert zero matrix to cell with no apical membrane
cell_no = 45;
mean_dist = 28.28; % distance of the cell from duct exit, can be found by running process_mesh_info.m .
dist = 100;
for i = 1:n_c
    cell_dist = abs(s_cell_prop{i}.mean_dist - mean_dist);
    % pick the index of the cell that is the closest to the skipped cell
    if cell_dist < dist
        min_i = i;
        dist = cell_dist;
    end
end

z_insert = zzz(:,min_i*n_c_remain+1:min_i*n_c_remain+n_c_remain);

insrt_idx = (cell_no - 1)*n_c_remain;
zzzz = [zzz(:,1:insrt_idx), z_insert, zzz(:,insrt_idx+1:end)];

%% get luminal water flow rates across time
% luminal water flow rate is not a direct output of ODE solver, hence we
% need to calculate it from the f_ODE_noMass function.

flowrate = zeros(n_t, n_l);
for i = 1:n_t
    [~, flowrate(i,:)] = f_ODE_noMass(t(i),z(i,:),P,s_cell_prop,s_lumen_prop,0,1,time_series);
end

%% subsampling with time, first 500s @ 0.1s step, second 100s at 1s step

time_sample = [1:5000,5001:10:10001];

t_sampled = t(time_sample);

zzzz = zzzz(time_sample,:);
save dynamic_data zzzz

flowrate = flowrate(time_sample,:);
save dynamic_flow flowrate

save lumen_prop s_lumen_prop
