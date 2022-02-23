function cell_prop = scale_parameters(cell_prop, P_i, P_s)
for i = 1:length(cell_prop)
    cell_struct = cell_prop{i};
    
    % scale the rates based on cell surface areas
    if cell_struct.type == 'S'
        P = P_s;
    elseif cell_struct.type == 'I'
        P = P_i;
    end
    
    A_A = cell_struct.api_area;
    A_B = cell_struct.baslat_area;
    n_c = cell_struct.n_c;
    
%     if cell_struct.mean_dist < 30
%         P.G_P_Cl = P.G_P_Cl*0.05;
%     end
    scaled_rates = parameter_scalling(P, A_A, A_B, n_c);
    
    scaled_rates.NHE.alpha_A = P.NHE.alpha_A;
    scaled_rates.NHE.alpha_B = P.NHE.alpha_B;
    scaled_rates.AE2.alpha_A = P.AE2.alpha_A;
    scaled_rates.AE2.alpha_B = P.AE2.alpha_B;
    scaled_rates.NBC.alpha_A = P.NBC.alpha_A;
    scaled_rates.NBC.alpha_B = P.NBC.alpha_B;
    
    cell_struct.scaled_rates = scaled_rates;
    cell_prop{i} = cell_struct;
end
end

function scaled_rates = parameter_scalling(P, A_A, A_B, n_c)
% the apical area used to scale the conductances G
A = 104.719755;
    
scaled_rates = struct;
scaled_rates.L_A    = n_c * P.L_A * A / A_A;
scaled_rates.L_B    = n_c * P.L_B * A / A_B;
scaled_rates.G_ENaC = n_c * P.G_ENaC * A / A_A;
scaled_rates.G_CFTR = n_c * P.G_CFTR * A / A_A;
scaled_rates.G_CaCC = n_c * P.G_CaCC * A / A_A;
scaled_rates.G_BK   = n_c * P.G_BK * A / A_A;
scaled_rates.G_K_B  = n_c * P.G_K_B * A / A_B;
scaled_rates.G_P_Na = n_c * P.G_P_Na * A / A_A;
scaled_rates.G_P_K  = n_c * P.G_P_K * A / A_A;
scaled_rates.G_P_Cl = n_c * P.G_P_Cl * A / A_A;
scaled_rates.NKA    = P.NKA;
scaled_rates.NKA.alpha_A = n_c * P.NKA.alpha_A * A / A_A;
scaled_rates.NKA.alpha_B = n_c * P.NKA.alpha_B * A / A_B;

% the protein molecules are multipled by the number of cells included in
% a simplified cell
scaled_rates.chi_C = n_c * P.chi_C;
end
