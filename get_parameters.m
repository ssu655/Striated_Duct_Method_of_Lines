function P = get_parameters(Conc,PSflow)

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

P = struct;
P.ConI = Conc.Int;
P.ConP = Conc.PS;

P.PSflow = PSflow;

% apical channels conductances 
P.G_ENaC = 2.5;

P.G_CFTR = 10;

P.G_BK = 6;

% basolateral channels conductances 
P.G_K_B = 0.5;

% apical or basolateral transporter rates
P.NBC = struct;
P.NBC.alpha = 100;
P.NBC.k5_p = -6e-1; % 1/s
P.NBC.k5_m = 1e8; % 1/s
P.NBC.k6_p = 1e8; % 1/s
P.NBC.k6_m = -1.9e-1; % 1/s

P.AE2 = struct;
P.AE2.alpha_A = 0.001;
P.AE2.alpha_B = 0.0001;
P.AE2.k3_p = 5.86; % 1/s
P.AE2.k3_m = 1.06e8; % 1/s
P.AE2.k4_p = 9.3e7; % 1/s
P.AE2.k4_m = 5.14; % 1/s

P.NHE = struct;
P.NHE.alpha_A = 0.0001;
P.NHE.alpha_B = 0.0001;
P.NHE.k1_p = 1.4e3; % 1/s
P.NHE.k1_m = 1.4e11; % 1/s
P.NHE.k2_p = 2.5e9; % 1/s
P.NHE.k2_m = 1.78e2; % 1/s

% CO2 permeability
P.p_CO = 50; % 1/s 

% CO2 bicarbonate buffering
P.buf = struct;
P.buf.k_p = 0.03; %/s
P.buf.k_m = 20; %/mMs

% sodium potassium pump rates
P.NKA = struct;
P.NKA.alpha_A = 0.7e-8; % mol/m2
P.NKA.alpha_B = 0.9e-8; % mol/m2

P.NKA.r = 1.305e-3; %mM-3s-1
P.NKA.beta = 0.647e-4; %mM-1

% paracellular conductances
P.G_P_Na = 0.1; %S/m2
P.G_P_K = 1; %S/m2
P.G_P_Cl = 1.5; %S/m2

% water permeability across membranes
P.L_A = 0;%0.6e1; % um/s
P.L_B = 0.6e1; % um/s

% universal physical constants
P.R = 8.13144621; % J/mol/K
P.T = 310; % K
P.F = 96485.3329; % C/mol
P.V = PSflow; % um^3 volumetric flow rate of primary saliva
P.V_w = 18e12; % um^3/mol partial molar mass of water

% osmolarity adjusting constants
P.chi_C = 4e-14; % mol (40 mM * 1000 um3  = xxx e-18 mol)
P.phi_A = 0.2; % mM (fong 2016)
P.phi_B = 10.92; % mM (Mangos 1972)
end