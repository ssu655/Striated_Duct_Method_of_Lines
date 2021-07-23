function P_list = get_parameters(Conc,PSflow,segment_meshes)

% P_list is a cell array of structures P: P_list = {P, P, ...} (1, n_s)
%
% The nested structure of P is as follows:
%
%    fields      sub-fields
% 
% P - ConI              (struct of Interstitium ion concentrations)
%   - ConP              (struct of primary saliva ion concentrations)
%   - PSflow
%   - G_ENaC
%   - G_CFTR
%   - G_K_B    
%   - G_BK      
%   - NBC       - alpha   
%               - k's
%   - AE2       - alpha_A
%               - alpha_B
%               - k's
%   - NHE       - alpha_A
%               - alpha_B
%               - k's
%   - p_CO
%   - buf       - k_m
%               - k_p
%   - NKA       - alpha
%               - r
%               - beta
%   - G_P_Na
%   - G_P_K
%   - G_P_Cl
%   - L_A
%   - L_B
%   - w_A
%   - chi_C
%   - phi_A
%   - phi_B
%   - R
%   - T
%   - F
%   - V
%   - A_L
%   - V_w
%

addpath('C:\Users\lingm\Dropbox\PhD\method_of_lines_multimesh\ini2struct')
INI = ini2struct('C:\Users\lingm\Dropbox\PhD\method_of_lines_multimesh\parms_default.ini');
param_s = INI.striated;
param_i = INI.intercalated;
param = INI.duct_common;

n_seg = length(segment_meshes.meshtypes);
P_list = cell(1,n_seg);

for i = 1:n_seg
    type = segment_meshes.meshtypes(i);
    
    P = struct;
    P.ConI = Conc.Int;
    P.ConP = Conc.PS;

    P.PSflow = PSflow;

    % apical channels conductances 
    if type == 1    % intercalated cell parameters
        P.G_ENaC = param_i.g_enac; 
    else            % striated cell parameters
        P.G_ENaC = param_s.g_enac;
    end

    P.G_CFTR = param.g_cftr;

    P.G_BK = param.g_bk;

    % basolateral channels conductances 
    P.G_K_B = param.g_k_b;

    % apical or basolateral transporter rates
    P.NBC = struct;
    P.NBC.alpha = param.nbc_alpha;
    P.NBC.k5_p = param.nbc_k5_p; % 1/s
    P.NBC.k5_m = param.nbc_k5_m; % 1/s
    P.NBC.k6_p = param.nbc_k6_p; % 1/s
    P.NBC.k6_m = param.nbc_k6_m; % 1/s

    P.AE2 = struct;
    P.AE2.alpha_A = param.ae2_alpha_a;
    P.AE2.alpha_B = param.ae2_alpha_b;
    P.AE2.k3_p = param.ae2_k3_p; % 1/s
    P.AE2.k3_m = param.ae2_k3_m; % 1/s
    P.AE2.k4_p = param.ae2_k4_p; % 1/s
    P.AE2.k4_m = param.ae2_k4_m; % 1/s

    P.NHE = struct;
    P.NHE.alpha_A = param.nhe_alpha_a;
    P.NHE.alpha_B = param.nhe_alpha_b;
    P.NHE.k1_p = param.nhe_k1_p; % 1/s
    P.NHE.k1_m = param.nhe_k1_m; % 1/s
    P.NHE.k2_p = param.nhe_k2_p; % 1/s
    P.NHE.k2_m = param.nhe_k2_m; % 1/s

    % CO2 permeability
    P.p_CO = param.p_co; % 1/s 

    % CO2 bicarbonate buffering
    P.buf = struct;
    P.buf.k_p = param.buf_k_p; %/s
    P.buf.k_m = param.buf_k_m; %/mMs

    % sodium potassium pump rates
    P.NKA = struct;
    if type == 1    % intercalated cell parameters
        P.NKA.alpha_A = param_i.nka_alpha_a; % mol/m2
        P.NKA.alpha_B = param_i.nka_alpha_b; % mol/m2
    else            % striated cell parameters
        P.NKA.alpha_A = param_s.nka_alpha_a; % mol/m2
        P.NKA.alpha_B = param_s.nka_alpha_b; % mol/m2
    end

    P.NKA.r = param.nka_r; %mM-3s-1
    P.NKA.beta = param.nka_beta; %mM-1

    % paracellular conductances
    if type == 1  
        P.G_P_Na = param_i.g_p_na; %S/m2
        P.G_P_K = param_i.g_p_k; %S/m2
        P.G_P_Cl = param_i.g_p_cl; %S/m2
    else
        P.G_P_Na = param_s.g_p_na; %S/m2
        P.G_P_K = param_s.g_p_k; %S/m2
        P.G_P_Cl = param_s.g_p_cl; %S/m2
    end
    
    % water permeability across membranes
    if type == 1  
        P.L_A = param_i.l_a; % um/s
        P.L_B = param_i.l_b; % um/s
    else
        P.L_A = param_s.l_a; % um/s
        P.L_B = param_s.l_b; % um/s
    end

    % universal physical constants
    P.R = 8.13144621; % J/mol/K
    P.T = 310; % K
    P.F = 96485.3329; % C/mol
    P.V = PSflow; % um^3 volumetric flow rate of primary saliva
    P.V_w = 18e12; % um^3/mol partial molar mass of water

    % osmolarity adjusting constants
    P.chi_C = param.chi_c; % mol (40 mM * 1000 um3  = xxx e-18 mol)
    P.phi_A = param.phi_a; % mM (fong 2016)
    P.phi_B = param.phi_b; % mM (Mangos 1972)
    
    % record the segment type in P
    P.seg_type = type;
    
    P_list{i} = P;
    
end
end
