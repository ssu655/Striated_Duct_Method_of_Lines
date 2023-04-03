function [P_i, P_s] = get_parameters(Conc, PSflow, parms_file)

% The nested structure of P is as follows:
%
%    fields      sub-fields
% 
% P - ConI              (struct of Interstitium ion concentrations)
%   - ConP              (struct of primary saliva ion concentrations)
%   - PSflow
%   - G_ENaC
%   - G_CFTR
%   - G_CaCC
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

addpath('ini2struct')
INI = ini2struct(parms_file);
param_s = INI.striated;
param_i = INI.intercalated;
param_c = INI.duct_common;

mergestructs = @(x,y) cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);

% striated parameters
param = mergestructs(param_s, param_c);
P_s = copy_param(param);

P_s.ConI = Conc.Int;
P_s.ConP = Conc.PS;
P_s.PSflow = PSflow;

% intercalated parameters
param = mergestructs(param_i, param_c);
P_i = copy_param(param);

P_i.ConI = Conc.Int;
P_i.ConP = Conc.PS;
P_i.PSflow = PSflow;

end

function P = copy_param(param)

a = 1;%40; % for parotid
% a = 60; % for SMG

% apical channels conductances 
P.G_ENaC = a*param.g_enac;% 2.5;

P.G_CFTR = a*param.g_cftr;% 10;

P.G_CaCC = a*param.g_cacc;% 10;

P.G_BK = a*param.g_bk;% 6;

% basolateral channels conductances 
P.G_K_B = a*param.g_k_b;% 0.5;
P.G_Na_B = a*param.g_na_b;% 0.5;
P.G_Cl_B = a*param.g_cl_b;% 0.5;

% apical or basolateral transporter rates
P.NBC = struct;
P.NBC.alpha_A = a*param.nbc_alpha_a; % 100;
P.NBC.alpha_B = a*param.nbc_alpha_b; % 100;
P.NBC.k5_p = param.nbc_k5_p; % -6e-1; % 1/s
P.NBC.k5_m = param.nbc_k5_m; % 1e8; % 1/s
P.NBC.k6_p = param.nbc_k6_p; % 1e8; % 1/s
P.NBC.k6_m = param.nbc_k6_m; % -1.9e-1; % 1/s

P.AE2 = struct;
P.AE2.alpha_A = a*param.ae2_alpha_a; % 0.001;
P.AE2.alpha_B = a*param.ae2_alpha_b; %0.0001;
P.AE2.k3_p = param.ae2_k3_p; %5.86; % 1/s
P.AE2.k3_m = param.ae2_k3_m; %1.06e8; % 1/s
P.AE2.k4_p = param.ae2_k4_p; %9.3e7; % 1/s
P.AE2.k4_m = param.ae2_k4_m; %5.14; % 1/s

P.NHE = struct;
P.NHE.alpha_A = a*param.nhe_alpha_a; %0.0001;
P.NHE.alpha_B = a*param.nhe_alpha_b; %0.0001;
P.NHE.k1_p = param.nhe_k1_p; %1.4e3; % 1/s
P.NHE.k1_m = param.nhe_k1_m; %1.4e11; % 1/s
P.NHE.k2_p = param.nhe_k2_p; %2.5e9; % 1/s
P.NHE.k2_m = param.nhe_k2_m; %1.78e2; % 1/s

% CO2 permeability
P.p_CO = a*param.p_co; %50; % 1/s 

% CO2 bicarbonate buffering
P.buf = struct;
P.buf.k_p = param.buf_k_p; %0.03; %/s
P.buf.k_m = param.buf_k_m; %20; %/mMs

% sodium potassium pump rates
P.NKA = struct;
P.NKA.alpha_A = a*param.nka_alpha_a; % 0.7e-8; % mol/m2
P.NKA.alpha_B = a*param.nka_alpha_b; % 0.9e-8; % mol/m2

P.NKA.r = param.nka_r; % 1.305e-3; %mM-3s-1
P.NKA.beta = param.nka_beta; % 0.647e-4; %mM-1

% paracellular conductances
P.G_P_Na = a*param.g_p_na; % 0.1; %S/m2
P.G_P_K = a*param.g_p_k; % 1; %S/m2
P.G_P_Cl = a*param.g_p_cl; % 1.5; %S/m2

% water permeability across membranes
P.L_A = param.l_a; % 0.6e1; % um/s
P.L_B = param.l_b; % 0.6e1; % um/s

% universal physical constants
P.R = 8.13144621; % J/mol/K
P.T = 310; % K
P.F = 96485.3329; % C/mol
P.V_w = 18e12; % um^3/mol partial molar mass of water

% osmolarity adjusting constants
P.chi_C = param.chi_c; % 4e-14; % mol (40 mM * 1000 um3  = xxx e-18 mol)
P.phi_A = param.phi_a; % 0.2; % mM (fong 2016)
P.phi_B = param.phi_b; % 10.92; % mM (Mangos 1972)
end