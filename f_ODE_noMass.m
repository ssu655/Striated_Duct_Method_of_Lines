function [dxdt,flow_rate] = f_ODE_noMass(t,x,P,cell_prop,lumen_prop,displ,dynamic,time_series)

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

n_c = length(cell_prop);  % number of cells
n_l = lumen_prop.n_disc;   % number of lumen segments

x_c = reshape(x(1 : n_c*9),9,[]);        % [9, n_c]
x_l = reshape(x(1 + n_c*9 : end),6,[]);  % [6, n_l]

if displ
    flux = initiate_flux(n_l, n_c, x_c);
end

% read the constant parameters from input struct P
if dynamic
    % primary saliva volume flow rate dependent on time
    if t<500
        Na_P = interp1(time_series.time,time_series.Na,t);
        K_P = interp1(time_series.time,time_series.K,t);
        Cl_P = interp1(time_series.time,time_series.Cl,t);
        HCO_P = interp1(time_series.time,time_series.HCO,t);
        H_P = interp1(time_series.time,time_series.H,t);
        pv = interp1(time_series.time,time_series.Q,t);
    else
        Na_P = P.ConP.Na;% mM
        K_P = P.ConP.K;
        Cl_P = P.ConP.Cl;
        HCO_P = P.ConP.HCO;
        H_P = P.ConP.H;
        pv = P.PSflow;
    end
else 
    % primary saliva volume flow rate independent on time
    Na_P = P.ConP.Na;% mM
    K_P = P.ConP.K;
    Cl_P = P.ConP.Cl;
    HCO_P = P.ConP.HCO;
    H_P = P.ConP.H;
    pv = P.PSflow;
end
% HCO_P = P.ConP.HCO;% 
% H_P = P.ConP.H;% 
CO_P = P.ConP.CO;% 

Na_B = P.ConI.Na;
K_B = P.ConI.K;
Cl_B = P.ConI.Cl; %
HCO_B = P.ConI.HCO; % 
H_B = P.ConI.H; % interstitium pH = 7.35
CO_B = P.ConI.CO;

R = P.R; % J/mol/K
T = P.T; % K
F = P.F;% C/mol

V_w = P.V_w; % um^3/mol 
A_L = lumen_prop.disc_X_area; %um^2 [1,n_l]
phi_A = P.phi_A; % mM
phi_B = P.phi_B; % mM

k1_p = P.NHE.k1_p; % 1/s
k1_m = P.NHE.k1_m; % 1/s
k2_p = P.NHE.k2_p; % 1/s
k2_m = P.NHE.k2_m; % 1/s

k3_p = P.AE2.k3_p; % 1/s
k3_m = P.AE2.k3_m; % 1/s
k4_p = P.AE2.k4_p; % 1/s
k4_m = P.AE2.k4_m; % 1/s

k5_p = P.NBC.k5_p; % 1/s
k5_m = P.NBC.k5_m; % 1/s
k6_p = P.NBC.k6_p; % 1/s
k6_m = P.NBC.k6_m; % 1/s

r_NKA = P.NKA.r; % mM-3s-1
beta_NKA = P.NKA.beta; % mM-1

p_CO = P.p_CO; % 1/s 
k_buf_p = P.buf.k_p; %/s
k_buf_m = P.buf.k_m; %/mMs

% setup a vector to record the rate of change of lumen fluid flow
dwAdt = zeros(1,n_l); 

% setup the ode rate of change matrices
dxcdt = zeros(size(x_c)); % [9,n_c]
dxldt = zeros(size(x_l)); % [6,n_l]

% loop through the cells to populate the rate of change for each cell/variable
for i = 1:n_c
    
    % first, read all the cell specific parameters
    cell_struct = cell_prop{i};
    
    G_ENaC = cell_struct.scaled_rates.G_ENaC;
    G_CFTR = cell_struct.scaled_rates.G_CFTR;
    G_CaCC = cell_struct.scaled_rates.G_CaCC;
    G_BK   = cell_struct.scaled_rates.G_BK;
    G_K_B  = cell_struct.scaled_rates.G_K_B;
    G_Cl_B = cell_struct.scaled_rates.G_Cl_B;
    G_Na_B = cell_struct.scaled_rates.G_Na_B;
    G_P_Na = cell_struct.scaled_rates.G_P_Na;
    G_P_K  = cell_struct.scaled_rates.G_P_K;
    G_P_Cl = cell_struct.scaled_rates.G_P_Cl;
    
    alpha_NKA_A = cell_struct.scaled_rates.NKA.alpha_A;
    alpha_NKA_B = cell_struct.scaled_rates.NKA.alpha_B;
    
    L_B = cell_struct.scaled_rates.L_B; % um/s 
    L_A = cell_struct.scaled_rates.L_A; % um/s 
    
    alpha_NHE_A = cell_struct.scaled_rates.NHE.alpha_A;
    alpha_NHE_B = cell_struct.scaled_rates.NHE.alpha_B;
    
    alpha_AE2_A = cell_struct.scaled_rates.AE2.alpha_A;
    alpha_AE2_B = cell_struct.scaled_rates.AE2.alpha_B;
    
    alpha_NBC_A = cell_struct.scaled_rates.NBC.alpha_A;
    alpha_NBC_B = cell_struct.scaled_rates.NBC.alpha_B;
    
    chi_C = cell_struct.scaled_rates.chi_C; % mol
    
    % read the cellular variables of this particular cell
    V_A   = x_c(1,i);
    V_B   = x_c(2,i);
    w_C   = x_c(3,i);
    Na_C  = x_c(4,i);
    K_C   = x_c(5,i);
    Cl_C  = x_c(6,i);
    HCO_C = x_c(7,i);
    H_C   = x_c(8,i);
    CO_C  = x_c(9,i);
    
    % read the lumenal variables of all lumen discs this cell interface with
    loc_disc = find(cell_struct.api_area_discs~=0);
    A_A = cell_struct.api_area;
    A_B = cell_struct.baslat_area;
    A_A_disc = cell_struct.api_area_discs(loc_disc); % um^2 [1,n_loc_disc]
    w_A = lumen_prop.disc_volume(loc_disc); % um^3 [1,n_loc_disc]
    
    Na_A  = x_l(1,loc_disc);
    K_A   = x_l(2,loc_disc);
    Cl_A  = x_l(3,loc_disc);
    HCO_A = x_l(4,loc_disc);
    H_A   = x_l(5,loc_disc);
    CO_A  = x_l(6,loc_disc);
    
    % water transport
    osm_c = chi_C./w_C*1e18; % osmolarity of cell due to proteins (chi)
    
    J_B = 1e-18*L_B.*V_w.*(Na_C + K_C + Cl_C + HCO_C + osm_c - Na_B - K_B - Cl_B - HCO_B - phi_B); % um/s 
    J_A = 1e-18*L_A.*V_w.*(Na_A + K_A + Cl_A + HCO_A + phi_A - Na_C - K_C - Cl_C - HCO_C - osm_c); % um/s [1, n_loc_disc]
    dwdt = A_B * J_B - sum(A_A_disc .* J_A); % um^3/s
    dwAdt(1,loc_disc) = dwAdt(1,loc_disc) + A_A_disc .* J_A; % um^3/s [1, n_loc_disc]
    
    % CDF C02 Diffusion 
    J_CDF_A = p_CO * (CO_C - CO_A).* w_C .* A_A_disc./A_A; % e-18 mol/s [1, n_loc_disc]
    J_CDF_B = p_CO * (CO_C - CO_B).* w_C; % e-18 mol/s
    
    % buf CO2 buffering
    J_buf_C = (k_buf_p*CO_C - k_buf_m.*HCO_C.*H_C).* w_C; % e-18 mol/s 
    J_buf_A =  k_buf_p*CO_A - k_buf_m.*HCO_A.*H_A; %mM/s [1,n_loc_disc]
    
    % NHE
    J_NHE_A = alpha_NHE_A*(k1_p*k2_p*Na_A.*H_C-k1_m*k2_m*Na_C.*H_A)./(k1_p*Na_A+k2_p*H_C+k1_m*H_A+k2_m*Na_C).*w_C.* A_A_disc./A_A; % e-18 mol/s [1,n_loc_disc]
    J_NHE_B = alpha_NHE_B*(k1_p*k2_p*Na_B *H_C-k1_m*k2_m*Na_C *H_B)./(k1_p*Na_B+k2_p*H_C+k1_m*H_B+k2_m*Na_C).*w_C; % e-18 mol/s 

    % AE2
    e_term = exp(F*0.001*V_A/(R*T));
    J_AE2_A = alpha_AE2_A*5e8.*(Cl_A.*HCO_C.^2 - e_term.*Cl_C.*HCO_A.^2)./(k3_p*Cl_A.^2+k4_p*HCO_C.^2+k3_m*HCO_A.^2+k4_m*Cl_C.^2).*w_C.* A_A_disc./A_A; % e-18 mol/s [1,n_loc_disc]
%     J_AE2_A = alpha_AE2_A*(k3_p*k4_p*Cl_A.*HCO_C - k3_m*k4_m*Cl_C.*HCO_A)./(k3_p*Cl_A+k4_p*HCO_C+k3_m*HCO_A+k4_m*Cl_C).*w_C.* A_A_disc./A_A; % e-18 mol/s [1,n_loc_disc]
    J_AE2_B = alpha_AE2_B*(k3_p*k4_p*Cl_B.*HCO_C - k3_m*k4_m*Cl_C.*HCO_B)./(k3_p*Cl_B+k4_p*HCO_C+k3_m*HCO_B+k4_m*Cl_C).*w_C; % e-18 mol/s 
    
    % NBC 
    J_NBC_B = alpha_NBC_B.*(k5_p*k6_p*Na_C.*HCO_C-k5_m*k6_m*Na_B*HCO_B)./(k5_p.*Na_C.*HCO_C+k6_p*k5_m+k6_m*Na_B*HCO_B).*w_C; % e-18 mol/s
    J_NBC_A = alpha_NBC_A.*(k5_p*k6_p*Na_C.*HCO_C-k5_m*k6_m.*Na_A.*HCO_A)./(k5_p.*Na_C.*HCO_C+k6_p*k5_m+k6_m.*Na_A.*HCO_A).*w_C.* A_A_disc./A_A; % e-18 mol/s [1,n_loc_disc]
    
    % CFTR
    V_A_Cl = 1e3*R*T/(-1*F).*log(Cl_A./Cl_C); % mV [1,n_loc_disc]
    I_CFTR = G_CFTR .* A_A_disc .* (V_A - V_A_Cl); % e-6 nA [1,n_loc_disc]
    I_CaCC = G_CaCC .* A_A_disc .* (V_A - V_A_Cl); % e-6 nA [1,n_loc_disc]
    
    % CFTR_B
    V_A_HCO = 1e3*R*T/((-1)*F).*log(HCO_A./HCO_C); % mV [1,n_loc_disc]
    I_CFTR_B = 0.25 * G_CFTR .* A_A_disc .* (V_A - V_A_HCO); % e-6 nA [1,n_loc_disc]
    
    % I_BK Apical
    V_A_K = 1e3*R*T/F.*log(K_A./K_C); % mV  [1,n_loc_disc]
    I_BK = G_BK .* A_A_disc .* (V_A - V_A_K); % e-6 nA 

    % I_Na_B Basolateral
    V_B_Na = 1e3*R*T/F.*log(Na_B./Na_C); % mV
    I_Na_B = G_Na_B * A_B .* (V_B - V_B_Na); % e-6 nA 
    
    % I_K_B Basolateral
    V_B_K = 1e3*R*T/F.*log(K_B./K_C); % mV
    I_K_B = G_K_B * A_B .* (V_B - V_B_K); % e-6 nA 
    
    % I_Cl_B Basolateral
    V_B_Cl = 1e3*R*T/(-1*F).*log(Cl_B./Cl_C); % mV
    I_Cl_B = G_Cl_B * A_B .* (V_B - V_B_Cl); % e-6 nA 

    % ENaC Apical
    V_A_Na = 1e3*R*T/F*log(Na_A./Na_C); % mV  [1,n_loc_disc]
    I_ENaC = G_ENaC .* A_A_disc .* (V_A - V_A_Na); % e-6 nA

    % NaKATPase, NKA 
    J_NKA_A = A_A_disc .* alpha_NKA_A * r_NKA .*(K_A.^2.*Na_C.^3)./(K_A.^2+beta_NKA*Na_C.^3); % 10^-12 mol/s [1,n_loc_disc]
    J_NKA_B = A_B .* alpha_NKA_B * r_NKA .*(K_B.^2.*Na_C.^3)./(K_B.^2+beta_NKA*Na_C.^3); % 10^-12 mol/s
    
    
    % Paracellular currents
    V_T      = V_A - V_B; % mV
    V_P_Na   = 1e3*R*T/F.*log(Na_A/Na_B); % mV [1,n_loc_disc]
    V_P_K    = 1e3*R*T/F.*log(K_A/K_B); % mV [1,n_loc_disc]
    V_P_Cl   = 1e3*R*T/(-F).*log(Cl_A/Cl_B); % mV [1,n_loc_disc]
    I_P_Na   = G_P_Na .* A_A_disc .* (V_T - V_P_Na); % e-6 nA [1,n_loc_disc]
    I_P_K    = G_P_K .* A_A_disc .* (V_T - V_P_K); % e-6 nA [1,n_loc_disc]
    I_P_Cl   = G_P_Cl .* A_A_disc .* (V_T - V_P_Cl); % e-6 nA [1,n_loc_disc]
    
    
    % V_A e-15 c/s
    dxcdt(1,i) = -100000.*(sum(F*J_NKA_A*1e3 + I_ENaC + I_BK + I_CFTR + I_CaCC + I_CFTR_B + I_P_Na + I_P_K + I_P_Cl - F.*J_AE2_A*1e-3));
    % V_B e-15 c/s
    dxcdt(2,i) = -100000.*(F*J_NKA_B*1e3 + I_K_B + I_Na_B + I_Cl_B - sum(I_P_Na + I_P_K + I_P_Cl));
    % w_C um^3
    dxcdt(3,i) = dwdt;
    % Na_C mM/s
    dxcdt(4,i) = -dwdt*Na_C/w_C + 1e3*(-sum(I_ENaC)./(F*w_C) - I_Na_B./(F*w_C)) - 1e6*(3*(J_NKA_B+sum(J_NKA_A))/w_C) + sum(J_NBC_A)/w_C + J_NBC_B/w_C + sum(J_NHE_A)/w_C + J_NHE_B/w_C;
    % K_C mM/s
    dxcdt(5,i) = -dwdt*K_C/w_C + 1e3*(-sum(I_BK)./(F*w_C) - I_K_B./(F*w_C)) + 1e6*(2*(J_NKA_B+sum(J_NKA_A))/w_C);
    % Cl_C mM/s
    dxcdt(6,i) = -dwdt*Cl_C/w_C + 1e3*(sum(I_CFTR)./(F*w_C) + sum(I_CaCC)./(F*w_C) + I_Cl_B./(F*w_C)) + sum(J_AE2_A)/w_C + J_AE2_B/w_C;
    % HCO_C mM/s
    dxcdt(7,i) = -dwdt*HCO_C/w_C + 1e3*(sum(I_CFTR_B)./(F*w_C)) + sum(J_NBC_A)/w_C + J_NBC_B/w_C - 2*sum(J_AE2_A)/w_C - J_AE2_B/w_C + J_buf_C/w_C;
    % H_C mM/s
    dxcdt(8,i) = -dwdt*H_C/w_C - sum(J_NHE_A)/w_C - J_NHE_B/w_C + J_buf_C/w_C;
    % CO_C mM/s
    dxcdt(9,i) = -dwdt*CO_C/w_C - sum(J_CDF_A)/w_C - J_CDF_B/w_C - J_buf_C/w_C;
    
    % Na_A mM/s
    dxldt(1,loc_disc) = dxldt(1,loc_disc) + 1e6*(3*J_NKA_A./w_A) + 1e3*(I_ENaC./(F*w_A)) + 1e3*(I_P_Na./(F*w_A)) - J_NHE_A./w_A - J_NBC_A./w_A;
    % K_A mM/s
    dxldt(2,loc_disc) = dxldt(2,loc_disc) - 1e6*(2*J_NKA_A./w_A) + 1e3*(I_BK./(F*w_A)) + 1e3*(I_P_K./(F*w_A));
    % Cl_A mM/s
    dxldt(3,loc_disc) = dxldt(3,loc_disc) + 1e3*(-I_CFTR./(F*w_A)) + 1e3*(-I_CaCC./(F*w_A)) + 1e3*(-I_P_Cl./(F*w_A)) - J_AE2_A./w_A;
    % HCO_A mM/s
    dxldt(4,loc_disc) = dxldt(4,loc_disc) + 1e3*(-I_CFTR_B./(F*w_A)) - J_NBC_A./w_A + 2*J_AE2_A./w_A + J_buf_A;
    % H_A mM/s
    dxldt(5,loc_disc) = dxldt(5,loc_disc) + J_NHE_A./w_A + J_buf_A;
    % CO_A mM/s
    dxldt(6,loc_disc) = dxldt(6,loc_disc) + J_CDF_A./w_A - J_buf_A;
    
    if displ
        flux.V_A_Na(loc_disc) = V_A_Na;
        flux.V_P_Na(loc_disc) = V_P_Na;
        flux.V_B_Na(i) = V_B_Na;
        flux.V_A_K(loc_disc) = V_A_K;
        flux.V_B_K(i) = V_B_K;
        flux.V_P_K(loc_disc) = V_P_K;
        flux.V_A_Cl(loc_disc) = V_A_Cl;
        flux.V_B_Cl(i) = V_B_Cl;
        flux.V_P_Cl(loc_disc) = V_P_Cl;
        flux.V_A_HCO(loc_disc) = V_A_HCO;
        flux.J_NHE_A(loc_disc) = flux.J_NHE_A(loc_disc) + J_NHE_A;
        flux.J_NHE_A_c(i) = sum(J_NHE_A);
        flux.J_NHE_B(i) = J_NHE_B;
        flux.J_AE2_A(loc_disc) = flux.J_AE2_A(loc_disc) + J_AE2_A;
        flux.J_AE2_A_c(i) = sum(J_AE2_A);
        flux.J_AE2_B(i) = J_AE2_B;
        flux.J_NBC_A(loc_disc) = flux.J_NBC_A(loc_disc) + J_NBC_A;
        flux.J_NBC_A_c(i) = sum(J_NBC_A);
        flux.J_NBC_B(i) = J_NBC_B;
        flux.J_NKA_A(loc_disc) = flux.J_NKA_A(loc_disc) + J_NKA_A;
        flux.J_NKA_A_c(i) = sum(J_NKA_A);
        flux.J_NKA_B(i) = J_NKA_B;
        flux.J_buf_A(loc_disc) = flux.J_buf_A(loc_disc) + J_buf_A.*w_A;
        flux.J_buf_A_c(i) = sum(J_buf_A.*w_A);
        flux.J_buf_C(i) = J_buf_C;
        flux.I_ENaC(loc_disc) = flux.I_ENaC(loc_disc) + I_ENaC;
        flux.I_ENaC_c(i) = sum(I_ENaC);
        flux.I_P_Na(loc_disc) = flux.I_P_Na(loc_disc) + I_P_Na;
        flux.I_P_Na_c(i) = sum(I_P_Na);
        flux.I_BK(loc_disc) = flux.I_BK(loc_disc) + I_BK;
        flux.I_BK_c(i) = sum(I_BK);
        flux.I_K_B(i) = I_K_B;
        flux.I_Cl_B(i) = I_Cl_B;
        flux.I_Na_B(i) = I_Na_B;
        flux.I_P_K(loc_disc) = flux.I_P_K(loc_disc) + I_P_K;
        flux.I_P_K_c(i) = sum(I_P_K);
        flux.I_CFTR(loc_disc) = flux.I_CFTR(loc_disc) + I_CFTR;
        flux.I_CFTR_c(i) = sum(I_CFTR);
        flux.I_CaCC(loc_disc) = flux.I_CaCC(loc_disc) + I_CaCC;
        flux.I_CaCC_c(i) = sum(I_CaCC);
        flux.I_P_Cl(loc_disc) = flux.I_P_Cl(loc_disc) + I_P_Cl;
        flux.I_P_Cl_c(i) = sum(I_P_Cl);
        flux.I_CFTR_B(loc_disc) = flux.I_CFTR_B(loc_disc) + I_CFTR_B;
        flux.I_CFTR_B_c(i) = sum(I_CFTR_B);
    end
end

% compute the fluid flow rate in the lumen
v_secreted = zeros(1,n_l); % um^3/s accumulated volume flow rates of secreted fluid
v_up = zeros(1,n_l); % um^3/s volume flow rate of fluid into each lumen disc
x_up = zeros(size(x_l)); 

disc_out = lumen_prop.disc_out_Vec;

for i=n_l:-1:1
    
    % find upstream disc(s) of disc i
    i_up = find(ismember(disc_out, i));
    
    % if no upstream disc, it is an acinus end disc, v_up = PSflow
    if isempty(i_up)
        v_up(i) = pv;
        v_secreted(i) = 0;
%         x_up(:,i) = cell2mat(struct2cell(P.ConP));
        x_up(:,i) = [Na_P; K_P; Cl_P; HCO_P; H_P; CO_P];
        
    else % i_up could be a vector due to i being a branching disc
        v_up(i) = sum(v_up(i_up));
        v_secreted(i) = sum(v_secreted(i_up) + dwAdt(i_up));
        x_up(:,i) = sum(x_l(:,i_up).*v_up(i_up)/sum(v_up(i_up)),2);
    end
end
v_up = v_up + v_secreted;
v = v_up + dwAdt;

% convert volumetric flow rate to linear flow speed
v = v./A_L; % um/s 
v_up = v_up./A_L; % um/s

% 1D finite difference discretisation of the lumen, backward differences scheme
for i = 1:6
    dxldt(i,:) = dxldt(i,:) + (v_up.*x_up(i,:) - v.*x_l(i,:))./lumen_prop.disc_length;
end

% flatten the matrix to a column vector
dxdt = [dxcdt(:); dxldt(:)];

flow_rate = v.*A_L;

% display for debugging and cross checking purposes
if displ
    fprintf('initial P.S. flow rate: %2.2f  um3 \n',(5*v_up(end)*A_L(end))) % um^3/s
    fprintf('final P.S. flow rate:   %2.2f  um3 \n',(v(1)*A_L(1))) % um^3/s
    fprintf('percentage:             %2.2f  ',(v(1)*A_L(1)-5*v_up(end)*A_L(end))/(5*v_up(end)*A_L(end))*100)
    
    IntPos = zeros(1,lumen_prop.n_disc);
    IntPos(1) = lumen_prop.disc_length(1);
    for i = 2:lumen_prop.n_disc
        out = lumen_prop.disc_out_Vec(i);
        IntPos(i) = lumen_prop.disc_length(i) + IntPos(out);
    end
    max_length = max(IntPos);
    IntPos = max_length - IntPos;
    % IntPos = IntPos(1:58);
    % y_l = y_l(:,1:58);

    CellPos = zeros(1,length(cell_prop));
    CellType = zeros(2, length(cell_prop));
    for i = 1:length(cell_prop)
        CellPos(i) = cell_prop{i}.mean_dist;
        if cell_prop{i}.type == "I"
            CellType(:,i) = [1,0];
        else
            CellType(:,i) = [0,1];
        end
    end
    CellPos = max_length - CellPos;
    
    figure(4)
    subplot(4,4,1)
    plot(CellPos,flux.V_A,'.')
    hold on
    plot(CellPos,flux.V_T,'.')
    plot(IntPos,flux.V_A_Na,'.')
    plot(CellPos,flux.V_B_Na,'.')
    plot(IntPos,flux.V_P_Na,'.')
    hold off
    ylabel('mV')
    legend('V_A','V_T','V_{A_{Na}}','V_{B_{Na}}','V_{P_{Na}}')
    
    subplot(4,4,5)
    plot(CellPos,flux.V_A,'.')
    hold on
    plot(CellPos,flux.V_T,'.')
    plot(CellPos,flux.V_B,'.')
    plot(IntPos,flux.V_A_K,'.')
    plot(CellPos,flux.V_B_K,'.')
    plot(IntPos,flux.V_P_K,'.')
    hold off
    ylabel('mV')
    legend('V_A','V_T','V_B','V_{A_K}','V_{B_K}','V_{P_K}')
    
    subplot(4,4,9)
    plot(CellPos,flux.V_A,'.')
    hold on
    plot(CellPos,flux.V_T,'.')
    plot(IntPos,flux.V_A_Cl,'.')
    plot(CellPos,flux.V_B_Cl,'.')
    plot(IntPos,flux.V_P_Cl,'.')
    hold off
    ylabel('mV')
    legend('V_A','V_T','V_{A_{Cl}}','V_{B_{Cl}}','V_{P_{Cl}}')
    
    subplot(4,4,13)
    plot(CellPos,flux.V_A,'.')
    hold on
    plot(IntPos,flux.V_A_HCO,'.')
    hold off
    ylabel('mV')
    legend('V_A','V_{A_{HCO}}')
    
    subplot(4,4,3)
    plot(IntPos,flux.J_NHE_A.*F.*1e-9,'.')
    hold on
    plot(CellPos,flux.J_NHE_B.*F.*1e-9,'.')
    hold off
    ylabel('Current nA')
    legend('J_{NHE_A}','J_{NHE_B}')
    
    subplot(4,4,7)
    plot(IntPos,flux.J_NKA_A*F*1e-3,'.')
    hold on
    plot(CellPos,flux.J_NKA_B*F*1e-3,'.')
    hold off
    ylabel('Current nA')
    legend('J_{NKA_A}','J_{NKA_B}')
    
    subplot(4,4,11)
    plot(IntPos,flux.J_AE2_A.*F.*1e-9,'.')
    hold on
    plot(CellPos,flux.J_AE2_B.*F.*1e-9,'.')
    hold off
    ylabel('Current nA')
    legend('J_{AE2_A}','J_{AE2_B}')
    
    subplot(4,4,15)
    plot(IntPos,flux.J_NBC_A.*F.*1e-9,'.')
    hold on
    plot(CellPos,flux.J_NBC_B.*F.*1e-9,'.')
    plot(IntPos,flux.J_buf_A.*F.*1e-9,'.')
    plot(CellPos,flux.J_buf_C.*F.*1e-9,'.')
    hold off
    ylabel('Current nA')
    legend('J_{NBC_A}','J_{NBC_B}','J_{buf_A}','J_{buf_C}') %'J_{CDF_A}','J_{CDF_B}'
    
    subplot(4,4,2)
    plot(IntPos,flux.I_ENaC*1e-6,'.')
    hold on
    plot(CellPos,flux.I_Na_B*1e-6,'.')
    plot(IntPos,flux.I_P_Na*1e-6,'.')
    hold off
    ylabel('Current nA')
    legend('I_{ENaC}','I_{Na_B}','I_{P_{Na}}')
    
    subplot(4,4,6)
    plot(IntPos,flux.I_BK*1e-6,'.')
    hold on
    plot(CellPos,flux.I_K_B*1e-6,'.')
    plot(IntPos,flux.I_P_K*1e-6,'.')
    hold off
    ylabel('Current nA')
    legend('I_{BK}','I_{K_B}','I_{P_K}')
    
    subplot(4,4,10)
    plot(IntPos,flux.I_CFTR*1e-6,'.')
    hold on
    plot(IntPos,flux.I_CaCC*1e-6,'.')
    plot(CellPos,flux.I_Cl_B*1e-6,'.')
    plot(IntPos,flux.I_P_Cl*1e-6,'.')
    hold off
    ylabel('Current nA')
    legend('I_{CFTR}','I_{CaCC}','I_{Cl_B}','I_{P_{Cl}}')
    
    subplot(4,4,14)
    plot(IntPos,flux.I_CFTR_B*1e-6,'.')
    ylabel('Current nA')
    legend('I_{CFTR_B}')
    
    subplot(4,4,4)
    plot(IntPos,flux.I_ENaC.*1e-6 - flux.J_NHE_A.*F.*1e-9 - flux.J_NBC_A.*F.*1e-9 + flux.I_P_Na.*1e-6 + 3*flux.J_NKA_A*1e-3*F,'.')
    legend('Na flux')
    subplot(4,4,8)
    plot(IntPos,flux.I_BK.*1e-6 + flux.I_P_K.*1e-6 - 2*flux.J_NKA_A*1e-3*F,'.')
    legend('K flux')
    subplot(4,4,12)
    plot(IntPos,flux.I_P_Cl.*1e-6 + flux.I_CFTR.*1e-6 + flux.J_AE2_A.*F*1e-9,'.')
    hold on
    plot(IntPos,flux.I_CFTR_B.*1e-6 - flux.J_AE2_A.*F.*1e-9 - flux.J_NBC_A.*F.*1e-9 - flux.J_buf_A.*F.*1e-9,'.')
    hold off
    legend('Cl flux','HCO flux')
    subplot(4,4,16)
    plot(IntPos,dwAdt,'.')
    legend('water flux')
   
%     c_idx = find(CellType(2,:));
%     c_idx = [1,12,23,34,35,43,45,51,54,56,60,75];
    c_idx = [44,45,49,42,71,64,79,13,68];
    t1 = 'WT:';

    % x axis plot range, proximal and distal
    lim_p = 53;
    lim_d = 145;
    figure(5)
    subplot(4,2,1)
    x = CellPos(c_idx);
    y = zeros(1,length(x));
    y(1,:) = flux.I_ENaC_c(c_idx).*1e-6;
    y(2,:) = flux.I_P_Na_c(c_idx).*1e-6;
    y(3,:) = -flux.J_NBC_A_c(c_idx).*F.*1e-9;
    y(4,:) = -flux.J_NHE_A_c(c_idx).*F.*1e-9;
    bar(x,y,2)
    % ylim([-0.3,0.003])
    xlim([lim_p,lim_d])
    title(strcat(t1,': Apical Na^+ into lumen '))
    ylabel('nA')
    legend('I_{ENaC}', 'I_{P_{Na}}', 'J_{NBC_A}', 'J_{NHE_A}')
    set(gca,'xtick',[])
    
    subplot(4,2,2)
    x = CellPos(c_idx);
    y = zeros(2,length(x));
    y(1,:) = -3.*flux.J_NKA_B(c_idx).*F*1e-3;
    y(2,:) = flux.J_NBC_B(c_idx).*F.*1e-9;
    y(3,:) = flux.J_NHE_B(c_idx).*F.*1e-9;
    y(4,:) = -flux.I_Na_B(c_idx).*1e-6;
    bar(x,y,2)
    % ylim([-0.4,0.1])
    xlim([lim_p,lim_d])
    title(strcat(t1,': Basolateral Na^+ into cell '))
    ylabel('nA')
    legend('J_{NKA_B}','J_{NBC_B}', 'J_{NHE_B}','I_{Na_B}')
    set(gca,'xtick',[])
    
    subplot(4,2,3)
    x = CellPos(c_idx);
    y = zeros(2,length(x));
    y(1,:) = flux.I_BK_c(c_idx).*1e-6;
    y(2,:) = flux.I_P_K_c(c_idx).*1e-6;
    bar(x,y,1.5)
    xlim([lim_p,lim_d])
    title(strcat(t1,': Apical K^+ into lumen '))
    ylabel('nA')
    legend('I_{BK}', 'I_{P_{K}}')
    set(gca,'xtick',[])

    subplot(4,2,4)
    x = CellPos(c_idx);
    y = zeros(2,length(x));
    y(1,:) = 2.*flux.J_NKA_B(c_idx).*F*1e-3;
    y(2,:) = -flux.I_K_B(c_idx).*1e-6;
    bar(x,y,1.5)
    xlim([lim_p,lim_d])
    title(strcat(t1,': Basolateral K^+ into cell '))
    ylabel('nA')
    legend('J_{NKA_B}','I_{K_B}')
    set(gca,'xtick',[])

    subplot(4,2,5)
    x = CellPos(c_idx);
    y = zeros(2,length(x));
    y(1,:) = -flux.I_CFTR_c(c_idx).*1e-6-flux.I_CaCC_c(c_idx).*1e-6;
    y(2,:) = -flux.I_P_Cl_c(c_idx).*1e-6;
    y(3,:) = -flux.J_AE2_A_c(c_idx).*F.*1e-9;
    bar(x,y,2)
    xlim([lim_p,lim_d])
    %ylim([-0.06,0.005])
    title(strcat(t1,': Apical Cl^- into lumen '))
    ylabel('nA') 
    legend('I_{CFTR}', 'I_{P_{Cl}}', 'J_{AE_A}')
    set(gca,'xtick',[])

    subplot(4,2,6)
    x = CellPos(c_idx);
    y = zeros(2,length(x));
    y(1,:) = flux.J_AE2_B(c_idx).*F.*1e-9;
    y(2,:) = flux.I_Cl_B(c_idx).*1e-6;
    bar(x,y,1.5)
    xlim([lim_p,lim_d])
    title(strcat(t1,': Basolateral Cl^- into cell ')) 
    ylabel('nA')
    legend('J_{AE_B}','I_{Cl_B}')
    set(gca,'xtick',[])

    subplot(4,2,7)
    x = CellPos(c_idx);
    y = zeros(2,length(x));
    y(1,:) = -flux.I_CFTR_B_c(c_idx).*1e-6;
    y(2,:) = flux.J_buf_A_c(c_idx).*F.*1e-9;
    y(3,:) = 2.*flux.J_AE2_A_c(c_idx).*F.*1e-9;
    y(4,:) = -flux.J_NBC_A_c(c_idx).*F.*1e-9;
    bar(x,y,2)
    xlim([lim_p,lim_d])
    % ylim([-0.04,0.04])
    title(strcat(t1,': Apical HCO_3^- into lumen '))
    ylabel('nA')
    xlabel('Duct entry                             Duct exit')
    legend('I_{CFTR_B}', 'J_{buf_A}', 'J_{AE_A}', 'J_{NBC_A}')
    set(gca,'xtick',[])

    subplot(4,2,8)
    x = CellPos(c_idx);
    y = zeros(2,length(x));
    y(1,:) = -flux.J_AE2_B(c_idx).*F.*1e-9;
    y(2,:) = -flux.J_buf_C(c_idx).*F.*1e-9;
    y(3,:) = flux.J_NBC_B(c_idx).*F.*1e-9;
    bar(x,y,2)
    xlim([lim_p,lim_d])
    title(strcat(t1,': Basolateral HCO_3^- into cell ')) 
    ylabel('nA')
    legend('J_{AE_B}','J_{buf_C}','J_{NBC_B}')
    xlabel('Duct entry                             Duct exit')
    set(gca,'xtick',[])

end
end

function flux = initiate_flux(n_l,n_c,x_c)
    flux = struct;
    flux.V_A = x_c(1,:);
    flux.V_B = x_c(2,:);
    flux.V_T = flux.V_A - flux.V_B;
    flux.V_A_Na = zeros(1, n_l);
    flux.V_B_Na = zeros(1, n_c);
    flux.V_P_Na = zeros(1, n_l);
    flux.V_A_K = zeros(1, n_l);
    flux.V_B_K = zeros(1, n_c);
    flux.V_P_K = zeros(1, n_l);
    flux.V_A_Cl = zeros(1, n_l);
    flux.V_B_Cl = zeros(1, n_c);
    flux.V_P_Cl = zeros(1, n_l);
    flux.V_A_HCO = zeros(1, n_l);
    flux.J_NHE_A = zeros(1, n_l);
    flux.J_NHE_A_c = zeros(1, n_c);
    flux.J_NHE_B = zeros(1, n_c);
    flux.J_AE2_A = zeros(1, n_l);
    flux.J_AE2_A_c = zeros(1, n_c);
    flux.J_AE2_B = zeros(1, n_c);
    flux.J_NBC_A = zeros(1, n_l);
    flux.J_NBC_A_c = zeros(1, n_c);
    flux.J_NBC_B = zeros(1, n_c);
    flux.J_NKA_A = zeros(1, n_l);
    flux.J_NKA_B = zeros(1, n_c);
    flux.I_ENaC = zeros(1, n_l);
    flux.I_ENaC_c = zeros(1, n_c);
    flux.I_P_Na = zeros(1, n_l);
    flux.I_P_Na_c = zeros(1, n_c);
    flux.I_BK = zeros(1, n_l);
    flux.I_BK_c = zeros(1, n_c);
    flux.I_K_B = zeros(1, n_c);
    flux.I_P_K = zeros(1, n_l);
    flux.I_P_K_c = zeros(1, n_c);
    flux.I_CFTR = zeros(1, n_l);
    flux.I_CFTR_c = zeros(1, n_c);
    flux.I_CaCC = zeros(1, n_l);
    flux.I_CaCC_c = zeros(1, n_c);
    flux.I_P_Cl = zeros(1, n_l);
    flux.I_P_Cl_c = zeros(1, n_c);
    flux.I_CFTR_B = zeros(1, n_l);
    flux.I_CFTR_B_c = zeros(1, n_c);
    flux.J_buf_A = zeros(1, n_l);
    flux.J_buf_A_c = zeros(1, n_c);
    flux.J_buf_C = zeros(1, n_c);
end