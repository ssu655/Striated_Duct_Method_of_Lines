function plot_fluxes(figure_no, CellPos, IntPos, flux, dwAdt)
    
    F = 96485.3329; % C/mol Faraday's constant 

    figure(figure_no)

    subplot(4,4,1)
    plot(CellPos,flux.V_A,'.')
    hold on
    plot(CellPos,flux.V_T,'.')
    plot(CellPos,flux.V_B,'.')
    plot(IntPos,flux.V_A_Na,'.')
    plot(CellPos,flux.V_B_Na,'.')
    plot(IntPos,flux.V_P_Na,'.')
    hold off
    ylabel('mV')
    legend('V_A','V_T','V_B','V_{A_{Na}}','V_{B_{Na}}','V_{P_{Na}}')
    
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
    plot(CellPos,flux.V_B,'.')
    plot(IntPos,flux.V_A_Cl,'.')
    plot(CellPos,flux.V_B_Cl,'.')
    plot(IntPos,flux.V_P_Cl,'.')
    hold off
    ylabel('mV')
    legend('V_A','V_T','V_B','V_{A_{Cl}}','V_{B_{Cl}}','V_{P_{Cl}}')
    
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
    plot(IntPos,flux.J_AE_A.*F.*1e-9,'.')
    hold on
    plot(CellPos,flux.J_AE_B.*F.*1e-9,'.')
    hold off
    ylabel('Current nA')
    legend('J_{AE_A}','J_{AE_B}')
    
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
    plot(IntPos,flux.I_P_Cl.*1e-6 + flux.I_CFTR.*1e-6 + flux.J_AE_A.*F*1e-9,'.')
    hold on
    plot(IntPos,flux.I_CFTR_B.*1e-6 - flux.J_AE_A.*F.*1e-9 - flux.J_NBC_A.*F.*1e-9 - flux.J_buf_A.*F.*1e-9,'.')
    hold off
    legend('Cl flux','HCO flux')
    
    subplot(4,4,16)
    plot(IntPos,dwAdt,'.')
    ylabel("\mum^3/s")
    legend('water flux')

    set(gcf,'position',[600,50,900,900])

%     figure(2)
%     subplot(3,1,1)
%     plot(IntPos,dwAdt,'.','MarkerSize',10)
%     ylabel("\mum^3/s")
%     title('Ductal Water Secretion')
%     set(gca,'xtick',[])
end