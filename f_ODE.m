function dxdt = f_ODE(t,x,P,cell_prop,lumen_prop,displ)

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

n_c = length(cell_prop);
n_l = lumen_prop.n_int;

x_c = reshape(x(1 : n_c*9),9,[]); % [9, n_c]
x_l = reshape(x(1 + n_c*9 : end),6,[]); % [6, n_l]

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
w_A = lumen_prop.volume;
L = lumen_prop.L;
A_L = lumen_prop.X_area;
chi_C = P.chi_C; % mol
phi_A = P.phi_A; % mol per lumen interval volumn 
phi_B = P.phi_B; % mM

alpha_NHE_A = P.NHE.alpha_A;
alpha_NHE_B = P.NHE.alpha_B;
k1_p = P.NHE.k1_p; % 1/s
k1_m = P.NHE.k1_m; % 1/s
k2_p = P.NHE.k2_p; % 1/s
k2_m = P.NHE.k2_m; % 1/s

alpha_AE2_A = P.AE2.alpha_A;
alpha_AE2_B = P.AE2.alpha_B;
k3_p = P.AE2.k3_p; % 1/s
k3_m = P.AE2.k3_m; % 1/s
k4_p = P.AE2.k4_p; % 1/s
k4_m = P.AE2.k4_m; % 1/s

alpha_NBC = P.NBC.alpha;
k5_p = P.NBC.k5_p; % 1/s
k5_m = P.NBC.k5_m; % 1/s
k6_p = P.NBC.k6_p; % 1/s
k6_m = P.NBC.k6_m; % 1/s

r_NKA = P.NKA.r; % mM-3s-1
beta_NKA = P.NKA.beta; % mM-1

p_CO = P.p_CO; % 1/s 
k_buf_p = P.buf.k_p; %/s
k_buf_m = P.buf.k_m; %/mMs

dwadt = zeros(1,n_l);

dxcdt = zeros(size(x_c)); % [9,n_c]
dxldt = zeros(size(x_l)); % [6,n_l]

for i = 1:n_c
    cell_struct = cell_prop{i};
    
    A_A = cell_struct.api_area;
    A_B = cell_struct.baslat_area;
    A_A_int = cell_struct.api_area_int;
    
    G_ENaC = cell_struct.scaled_rates.G_ENaC;
    G_CFTR = cell_struct.scaled_rates.G_CFTR;
    G_BK   = cell_struct.scaled_rates.G_BK;
    G_K_B  = cell_struct.scaled_rates.G_K_B;
    G_P_Na = cell_struct.scaled_rates.G_P_Na;
    G_P_K  = cell_struct.scaled_rates.G_P_K;
    G_P_Cl = cell_struct.scaled_rates.G_P_Cl;
    
    alpha_NKA_A = cell_struct.scaled_rates.NKA.alpha_A;
    alpha_NKA_B = cell_struct.scaled_rates.NKA.alpha_B;
    
    V_A   = x_c(1,i);
    V_B   = x_c(2,i);
    w_C   = x_c(3,i);
    Na_C  = x_c(4,i);
    K_C   = x_c(5,i);
    Cl_C  = x_c(6,i);
    HCO_C = x_c(7,i);
    H_C   = x_c(8,i);
    CO_C  = x_c(9,i);
    
    loc_int = cell_struct.loc_int;
    
    Na_A  = x_l(1,loc_int);
    K_A   = x_l(2,loc_int);
    Cl_A  = x_l(3,loc_int);
    HCO_A = x_l(4,loc_int);
    H_A   = x_l(5,loc_int);
    CO_A  = x_l(6,loc_int);
    
    % water transport
    osm_c = chi_C./w_C*1e18; % osmolarity of cell due to proteins (chi)
    L_B = cell_struct.scaled_rates.L_B; % um/s 
    L_A = cell_struct.scaled_rates.L_A; % um/s 
    
    J_B = 1e-18*L_B.*V_w.*(Na_C + K_C + Cl_C + HCO_C + osm_c - Na_B - K_B - Cl_B - HCO_B - phi_B); % um/s 
    J_A = 1e-18*L_A.*V_w.*(Na_A + K_A + Cl_A + HCO_A + phi_A - Na_C - K_C - Cl_C - HCO_C - osm_c); % um/s [1, n_loc_int]
    dwdt = A_B * J_B - sum(A_A_int .* J_A); % um^3/s
    dwadt(1,loc_int) = dwadt(1,loc_int) + A_A_int .* J_A; % um^3/s [1, n_loc_int]
    
    % CDF C02 Diffusion 
    J_CDF_A = p_CO * (CO_C - CO_A).* w_C .* A_A_int./A_A; % e-18 mol/s [1, n_loc_int]
    J_CDF_B = p_CO * (CO_C - CO_B).* w_C; % e-18 mol/s
    
    % buf CO2 buffering
    J_buf_C = (k_buf_p*CO_C - k_buf_m.*HCO_C.*H_C).* w_C; % e-18 mol/s 
    J_buf_A =  k_buf_p*CO_A - k_buf_m.*HCO_A.*H_A; %mM/s [1,n_loc_int]
    
    % NHE
    J_NHE_A = alpha_NHE_A*(k1_p*k2_p*Na_A.*H_C-k1_m*k2_m*Na_C.*H_A)./(k1_p*Na_A+k2_p*H_C+k1_m*H_A+k2_m*Na_C).*w_C.* A_A_int./A_A; % e-18 mol/s [1,n_loc_int]
    J_NHE_B = alpha_NHE_B*(k1_p*k2_p*Na_B *H_C-k1_m*k2_m*Na_C *H_B)./(k1_p*Na_B+k2_p*H_C+k1_m*H_B+k2_m*Na_C).*w_C; % e-18 mol/s 

    % AE2
    J_AE2_A = alpha_AE2_A*(k3_p*k4_p*Cl_A.*HCO_C - k3_m*k4_m*Cl_C.*HCO_A)./(k3_p*Cl_A+k4_p*HCO_C+k3_m*HCO_A+k4_m*Cl_C).*w_C.* A_A_int./A_A; % e-18 mol/s [1,n_loc_int]
    J_AE2_B = alpha_AE2_B*(k3_p*k4_p*Cl_B.*HCO_C - k3_m*k4_m*Cl_C.*HCO_B)./(k3_p*Cl_B+k4_p*HCO_C+k3_m*HCO_B+k4_m*Cl_C).*w_C; % e-18 mol/s 

    % NBC Basolateral
    J_NBC = alpha_NBC.*(k5_p*k6_p*Na_C.*HCO_C-k5_m*k6_m*Na_B*HCO_B)./(k5_p.*Na_C.*HCO_C+k6_p*k5_m+k6_m*Na_B*HCO_B).*w_C; % e-18 mol/s

    % CFTR Apical 
    V_A_Cl = 1e3*R*T/(-1*F).*log(Cl_A./Cl_C); % mV [1,n_loc_int]
    I_CFTR = G_CFTR .* A_A_int .* (V_A - V_A_Cl); % e-6 nA [1,n_loc_int]
    
    % CFTR_B Apical
    V_A_HCO = 1e3*R*T/((-1)*F).*log(HCO_A./HCO_C); % mV [1,n_loc_int]
    I_CFTR_B = 0.25 * G_CFTR .* A_A_int .* (V_A - V_A_HCO); % e-6 nA [1,n_loc_int]
    
    % I_BK Apical
    V_A_K = 1e3*R*T/F.*log(K_A./K_C); % mV  [1,n_loc_int]
    I_BK = G_BK .* A_A_int .* (V_A - V_A_K); % e-6 nA 

    % I_K_B Basolateral
    V_B_K = 1e3*R*T/F.*log(K_B./K_C); % mV
    I_K_B = G_K_B * A_B .* (V_B - V_B_K); % e-6 nA 

    % ENaC Apical
    V_A_Na = 1e3*R*T/F*log(Na_A./Na_C); % mV  [1,n_loc_int]
    I_ENaC = G_ENaC .* A_A_int .* (V_A - V_A_Na); % e-6 nA

    % NaKATPase, NKA 
    J_NKA_A = A_A_int .* alpha_NKA_A * r_NKA .*(K_A.^2.*Na_C.^3)./(K_A.^2+beta_NKA*Na_C.^3); % 10^-12 mol/s [1,n_loc_int]
    J_NKA_B = A_B .* alpha_NKA_B * r_NKA .*(K_B.^2.*Na_C.^3)./(K_B.^2+beta_NKA*Na_C.^3); % 10^-12 mol/s
    
    
    % Paracellular currents
    V_T      = V_A - V_B; % mV
    V_P_Na   = 1e3*R*T/F.*log(Na_A/Na_B); % mV [1,n_loc_int]
    V_P_K    = 1e3*R*T/F.*log(K_A/K_B); % mV [1,n_loc_int]
    V_P_Cl   = 1e3*R*T/(-F).*log(Cl_A/Cl_B); % mV [1,n_loc_int]
    I_P_Na   = G_P_Na .* A_A_int .* (V_T - V_P_Na); % e-6 nA [1,n_loc_int]
    I_P_K    = G_P_K .* A_A_int .* (V_T - V_P_K); % e-6 nA [1,n_loc_int]
    I_P_Cl   = G_P_Cl .* A_A_int .* (V_T - V_P_Cl); % e-6 nA [1,n_loc_int]
    
   
    % V_A e-15 c/s
    dxcdt(1,i) = sum(F*J_NKA_A*1e3 + I_ENaC + I_BK + I_CFTR + I_CFTR_B + I_P_Na + I_P_K + I_P_Cl);
    % V_B e-15 c/s
    dxcdt(2,i) = F*J_NKA_B*1e3 + I_K_B - sum(I_P_Na + I_P_K + I_P_Cl);
    % w_C um^3
    dxcdt(3,i) = dwdt;
    % Na_C e-18 mol/s
    dxcdt(4,i) = + 1e3*(-sum(I_ENaC)./(F)) - 1e6*(3*(J_NKA_B+sum(J_NKA_A))) + J_NBC + sum(J_NHE_A) + J_NHE_B;
    % K_C e-18 mol/s
    dxcdt(5,i) = + 1e3*(-sum(I_BK)./(F) - I_K_B./(F)) + 1e6*(2*(J_NKA_B+sum(J_NKA_A)));
    % Cl_C e-18 mol/s
    dxcdt(6,i) = + 1e3*(sum(I_CFTR)./(F)) + sum(J_AE2_A) + J_AE2_B;
    % HCO_C e-18 mol/s
    dxcdt(7,i) = + 1e3*(sum(I_CFTR_B)./(F)) + J_NBC - sum(J_AE2_A) - J_AE2_B + J_buf_C;
    % H_C e-18 mol/s
    dxcdt(8,i) = - sum(J_NHE_A) - J_NHE_B + J_buf_C;
    % CO_C e-18 mol/s
    dxcdt(9,i) = - sum(J_CDF_A) - J_CDF_B - J_buf_C;
    
    % Na_A mM/s
    dxldt(1,loc_int) = dxldt(1,loc_int) + 1e6*(3*J_NKA_A./w_A) + 1e3*(I_ENaC./(F*w_A)) + 1e3*(I_P_Na./(F*w_A)) - J_NHE_A./w_A;
    % K_A mM/s
    dxldt(2,loc_int) = dxldt(2,loc_int) - 1e6*(2*J_NKA_A./w_A) + 1e3*(I_BK./(F*w_A)) + 1e3*(I_P_K./(F*w_A));
    % Cl_A mM/s
    dxldt(3,loc_int) = dxldt(3,loc_int) + 1e3*(-I_CFTR./(F*w_A)) + 1e3*(-I_P_Cl./(F*w_A)) - J_AE2_A./w_A;
    % HCO_A mM/s
    dxldt(4,loc_int) = dxldt(4,loc_int) + 1e3*(-I_CFTR_B./(F*w_A)) + J_AE2_A./w_A + J_buf_A;
    % H_A mM/s
    dxldt(5,loc_int) = dxldt(5,loc_int) + J_NHE_A./w_A + J_buf_A;
    % CO_A mM/s
    dxldt(6,loc_int) = dxldt(6,loc_int) + J_CDF_A./w_A - J_buf_A;
end

v = ones(1,n_l) * P.PSflow; % um^3/s volume flow rate of fluid out of each lumen compartment
v_up = ones(1,n_l) * P.PSflow; % um^3/s volume flow rate of fluid into each lumen compartment

x_up = zeros(size(x_l)); % construct a matrix for the upstream apical concentration for each cell
x_up(:,1) = cell2mat(struct2cell(P.ConP));

for i=1:n_l
    v(i) = v(i) + sum(dwadt(1,1:i));
end

if n_l>1 % if there are more than one intervals
    v_up(2:n_l) = v(1:n_l-1);
    x_up(:,2:n_l) = x_l(:,1:n_l-1);
end

v = v./A_L; % um/s convert volume flow rate to linear flow speed
v_up = v_up./A_L; % um/s convert volume flow rate to linear flow speed

for i = 1:6
    dxldt(i,:) = dxldt(i,:) + (v_up.*x_up(i,:) - v.*x_l(i,:))./L;
end

dxdt = [dxcdt(:); dxldt(:)];
% 
% if ~displ
%     hold on
%     plot(t,I_ENaC*1e-6,'.b')
%     plot(t,I_BK*1e-6,'.r')
%     plot(t,I_K_B*1e-6,'.y')
%     plot(t,I_CFTR*1e-6,'.c')
%     plot(t,J_NKA_A*F*1e-3,'.k')
%     plot(t,J_NKA_B*F*1e-3,'.m')
%     plot(t,J_AE2_A*F*w_C*1e-9,'+b')
%     plot(t,J_AE2_B*F*w_C*1e-9,'+r')
%     % plot(t,J_hyd_C*F*x(3)*1e-9,'*r')
%     plot(t,J_buf_C*F*w_C*1e-9,'*b')
%     plot(t,J_NHE_A*F*w_C*1e-9,'or')
%     plot(t,J_NHE_B*F*w_C*1e-9,'ok')
% end

if displ
    fprintf('initial P.S. flow rate: %2.2f  um3 \n',(v_up(1)*A_L)) % um^3/s
    fprintf('final P.S. flow rate:   %2.2f  um3 \n',(v(end)*A_L)) % um^3/s
    fprintf('percentage:             %2.2f  ',(v(end)-v_up(1))/v_up(1)*100)
    fprintf('\n')
    fprintf('I_ENaC:       %.8d nA  \n',I_ENaC*1e-6)
    fprintf('I_BK:         %.8d nA  \n',I_BK*1e-6)
    fprintf('I_K_B:        %.8d nA  \n',I_K_B*1e-6)
    fprintf('J_NKA_A:      %.8d nA  \n',J_NKA_A*F*1e-3)
    fprintf('J_NKA_B:      %.8d nA  \n',J_NKA_B*F*1e-3)
    fprintf('I_CFTR:       %.8d nA  \n',I_CFTR*1e-6)
    fprintf('I_CFTR_B:     %.8d nA  \n',I_CFTR_B*1e-6)
    fprintf('J_AE2_A:      %.8d nA  \n',J_AE2_A.*F.*1e-9)
    fprintf('J_AE2_B:      %.8d nA  \n',J_AE2_B.*F.*1e-9)
    fprintf('J_NBC:        %.8d nA  \n',J_NBC.*F.*1e-9)
    fprintf('J_NHE_A:      %.8d nA  \n',J_NHE_A.*F.*1e-9)
    fprintf('J_NHE_B:      %.8d nA  \n',J_NHE_B.*F.*1e-9)
    fprintf('J_CDF_A:      %.8d nA  \n',J_CDF_A.*F.*1e-9)
    fprintf('J_CDF_B:      %.8d nA  \n',J_CDF_B.*F.*1e-9)
    fprintf('J_buf_A:      %.8d nA  \n',J_buf_A.*F.*w_A.*1e-9)
    fprintf('J_buf_C:      %.8d nA  \n',J_buf_C.*F.*1e-9)
    fprintf('J_A:          %.8d nA  \n',J_A)
    fprintf('J_B:          %.8d nA  \n',J_B)
    fprintf(' \n')
    fprintf('V_A_K:     %.8d mV \n', V_A_K)
    fprintf('V_B_K:     %.8d mV \n', V_B_K)
    fprintf('V_A_Cl:    %.8d mV \n', V_A_Cl)
    fprintf('V_A_Na:    %.8d mV \n', V_A_Na)
    fprintf('V_A:       %.8d mV \n', x_c(1,:))
    fprintf('V_B:       %.8d mV \n', x_c(2,:))
    fprintf('V_T:       %.8d mV \n', V_T)
    fprintf(' \n')
    fprintf('V_P_Na:    %.8d mV \n', V_P_Na)
    fprintf('V_P_K:     %.8d mV \n', V_P_K)
    fprintf('V_P_Cl:    %.8d mV \n', V_P_Cl)
    fprintf('I_P_Na:    %.8d nA \n',I_P_Na*1e-6)
    fprintf('I_P_K:     %.8d nA \n',I_P_K*1e-6)
    fprintf('I_P_Cl:    %.8d nA \n',I_P_Cl*1e-6)   
    fprintf(' \n') 
    fprintf('Na flux A: %.8d nA \n',I_ENaC*1e-6 - J_NHE_A*F*1e-9 + I_P_Na*1e-6 + 3*J_NKA_A*1e-3*F)
    fprintf('K flux A:  %.8d nA \n',I_BK*1e-6 + I_P_K*1e-6 - 2*J_NKA_A*1e-3*F)
    fprintf('Cl flux A: %.8d nA \n', I_P_Cl*1e-6 + I_CFTR*1e-6 + J_AE2_A*F*1e-9)
    fprintf('HC flux A: %.8d nA \n', I_CFTR_B*1e-6 - J_AE2_A*F*1e-9 - J_buf_A.*F.*w_A.*1e-9)
    
    
%     CellPos = zeros(1,length(cell_prop));
%     for i = 1:length(cell_prop)
%         CellPos(i) = cell_prop{i}.centroid(3);
%     end
%     [CellPos,I] = sort(CellPos);
%     IntPos = lumen_prop.segment(1:end-1);
% 
%     
%     figure
%     subplot(4,4,1)
%     plot(IntPos,V_A_Na)
%     hold on
%     plot(IntPos,V_P_Na)
%     plot(IntPos,V_A_f,'--')
%     plot(IntPos,V_T_f,'--')
%     hold off
%     ylabel('mV')
%     legend('V_{A_{Na}}','V_{P_{Na}}','V_A','V_T')
%     
%     subplot(4,4,5)
%     plot(IntPos,V_A_K)
%     hold on
%     plot(CellPos,V_B_K(I))
%     plot(IntPos,V_P_K)
%     plot(IntPos,V_A_f,'--')
%     plot(CellPos,V_B(I),'--')
%     plot(IntPos,V_T_f,'--')
%     hold off
%     ylabel('mV')
%     legend('V_{A_K}','V_{B_K}','V_{P_K}','V_A','V_B','V_T')
%     
%     subplot(4,4,9)
%     plot(IntPos,V_A_Cl)
%     hold on
%     plot(IntPos,V_P_Cl)
%     plot(IntPos,V_A_f,'--')
%     plot(IntPos,V_T_f,'--')
%     hold off
%     ylabel('mV')
%     legend('V_{A_{Cl}}','V_{P_{Cl}}','V_A','V_T')
%     
%     subplot(4,4,13)
%     plot(IntPos,V_A_HCO)
%     hold on
%     plot(IntPos,V_A_f,'--')
%     hold off
%     ylabel('mV')
%     legend('V_{A_{HCO}}','V_A')
%     
%     subplot(4,4,3)
%     plot(CellPos,J_NHE_A_c.*F.*1e-9)
%     hold on
%     plot(CellPos,J_NHE_B.*F.*1e-9)
%     hold off
%     ylabel('Current nA')
%     legend('J_{NHE_A}','J_{NHE_B}')
%     
%     subplot(4,4,7)
%     plot(CellPos,J_NKA_A_c*F*1e-3)
%     hold on
%     plot(CellPos,J_NKA_B*F*1e-3)
%     hold off
%     ylabel('Current nA')
%     legend('J_{NKA_A}','J_{NKA_B}')
%     
%     subplot(4,4,11)
%     plot(CellPos,J_AE2_A_c.*F.*1e-9)
%     hold on
%     plot(CellPos,J_AE2_B.*F.*1e-9)
%     hold off
%     ylabel('Current nA')
%     legend('J_{AE2_A}','J_{AE2_B}')
%     
%     subplot(4,4,15)
%     plot(CellPos,J_NBC.*F.*1e-9)
%     hold on
%     plot(CellPos,J_buf_A.*F.*w_A.*1e-9*ConMat')
%     plot(CellPos,J_buf_C.*F.*1e-9)
%     hold off
%     ylabel('Current nA')
%     legend('J_{NBC}','J_{buf_A}','J_{buf_C}') %'J_{CDF_A}','J_{CDF_B}'
%     
%     subplot(4,4,2)
%     plot(CellPos,I_ENaC_c*1e-6)
%     hold on
%     plot(CellPos,I_P_Na_c*1e-6,'--')
%     hold off
%     ylabel('Current nA')
%     legend('I_{ENaC}','I_{P_{Na}}')
%     
%     subplot(4,4,6)
%     plot(CellPos,I_BK_c*1e-6)
%     hold on
%     plot(CellPos,I_K_B*1e-6)
%     plot(CellPos,I_P_K_c*1e-6,'--')
%     hold off
%     ylabel('Current nA')
%     legend('I_{BK}','I_{K_B}','I_{P_K}')
%     
%     subplot(4,4,10)
%     plot(CellPos,I_CFTR_c*1e-6)
%     hold on
%     plot(CellPos,I_P_Cl_c*1e-6,'--')
%     hold off
%     ylabel('Current nA')
%     legend('I_{CFTR}','I_{P_{Cl}}')
%     
%     subplot(4,4,14)
%     plot(CellPos,I_CFTR_B_c*1e-6)
%     ylabel('Current nA')
%     legend('I_{CFTR_B}')
%     
%     subplot(4,4,4)
%     plot(CellPos,I_ENaC_c.*1e-6 - J_NHE_A_c.*F.*1e-9 + I_P_Na_c.*1e-6 + 3*J_NKA_A_c*1e-3*F)
%     legend('Na flux')
%     subplot(4,4,8)
%     plot(CellPos,I_BK_c.*1e-6 + I_P_K_c.*1e-6 - 2*J_NKA_A_c*1e-3*F)
%     legend('K flux')
%     subplot(4,4,12)
%     plot(CellPos,I_P_Cl_c.*1e-6 + I_CFTR_c.*1e-6 + J_AE2_A_c.*F*1e-9)
%     legend('Cl flux')
%     subplot(4,4,16)
%     plot(CellPos,I_CFTR_B_c.*1e-6 - J_AE2_A_c.*F.*1e-9 - J_buf_A.*F.*w_A.*1e-9*ConMat')
%     legend('HCO flux')
end
end