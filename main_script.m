%%
% Author: Shan Su
% Date: July 12th 2021

% This is the main script for running and setting up the parotid gland
% duct fluid transport model using some mesh files as generated
% from experimental scan of mouse parotid gland tissue. The cell geometry
% and lumen geometry are all extracted from the ply mesh files.

% In this script, user can specify the volumetric flow rates of the 
% primary saliva and the discretisation interval length of the lumenal 
% compartment in um. 

% clear

%% Model input setup

L_int = 2; % um length of lumen discretisation interval
PSflow = 150/10; % um3/s volumetric primary saliva flow rate

fields = {'Na'; 'K'; 'Cl'; 'HCO'; 'H'; 'CO'};
Int = [140.2; 5.3; 102.6; 24.7; 1000*10^(-7.35); 1.28]; % concentration of interstitium
PS = [143.5; 5.2; 114.5; 34.2; 1000*10^(-7.35); 1.28];  % concentration of Primary Saliva
CIC = [17; 140; 22; 75; 1000*10^(-7.35); 1.28];  % cellular initial concentration
LIC = [143.5; 5.2; 114.5; 34.2; 1000*10^(-7.35); 1.28]; % lumenal initial concentration
Conc = struct;
Conc.Int = cell2struct(num2cell(Int),fields);
Conc.PS = cell2struct(num2cell(PS),fields);
Conc.CIC = cell2struct(num2cell(CIC),fields);
Conc.LIC = cell2struct(num2cell(LIC),fields);

%% Collect mesh files
segment_meshes = gather_mesh_files();

%% Parameter structure setup

P_list = get_parameters(Conc,PSflow,segment_meshes);

%% Get Mesh Info

seg_prop = process_meshes(P_list, L_int, segment_meshes);

%% Initial condition setup

x = setup_IC(Conc, seg_prop);

% f_ODE(1,x,P,0)

%% Setup and solve the ODE

% f_ODE(1,x,P,cell_prop,lumen_prop,1);

tspan = [0,500000];

% ===========================================
% % run the version of ODE with mass matrix
% tic
% options = odeset('Mass', @(t,x) mass(t,x,cell_prop,lumen_prop), 'MassSingular', 'yes');
% [t,y] = ode15s(@(t,y) f_ODE(t,y,P,cell_prop,lumen_prop,0), tspan, x, options);
% toc

% ==========================================
% run the version of ODE without mass matrix
tic
[t,y] = ode15s(@(t,y) f_ODE_noMass(t,y,P_list,seg_prop,0), tspan, x);
toc

% f_ODE(1,x,P,cell_prop,lumen_prop,1);

%% Sanity checks

[x_c, x_l] = reshape_variables(x, seg_prop);
[y_c, y_l] = reshape_variables(y(end,:), seg_prop);
% n_c = length(cell_prop);
% x_c = reshape(x(1 : n_c*9),9,[]); %[9, n_c]
% x_l = reshape(x(1+n_c*9 : end),6,[]); %[6, n_l]
% y_c = reshape(y(end,1 : n_c*9),9,[]); %[9, n_c]
% y_l = reshape(y(end,1+n_c*9 : end),6,[]); %[6, n_l]

% disp('Electroneutrality check cell: initial (mol)')
% disp((x_c(4,:)+x_c(5,:)-x_c(6,:)-x_c(7,:)+x_c(8,:)).*x_c(3,:)*1e-18-1.5*P.chi_C)
% 
% disp('Electroneutrality check cell: initial (mM)')
% disp((x_c(4,:)+x_c(5,:)-x_c(6,:)-x_c(7,:)+x_c(8,:)-1.5*P.chi_C./x_c(3,:)*1e18))
% 
% disp('Electroneutrality check cell: final (mol)')
% disp((y_c(4,:)+y_c(5,:)-y_c(6,:)-y_c(7,:)+y_c(8,:)).*y_c(3,:)*1e-18-1.5*P.chi_C)
% 
% disp('Electroneutrality check cell: final (mM)')
% disp((y_c(4,:)+y_c(5,:)-y_c(6,:)-y_c(7,:)+y_c(8,:)-1.5*P.chi_C./y_c(3,:)*1e18))
% 
% disp('Electroneutrality check lumen: initial (mM)')
% disp((x_l(2,1)+x_l(3,1)-x_l(4,1)-x_l(5,1)+x_l(6,1)).*lumen_prop.volume*1e-18)
% disp('Electroneutrality check lumen: final (mM)')
% disp((y_l(2,:)+y_l(3,:)-y_l(4,:)-y_l(5,:)+y_l(6,:)).*lumen_prop.volume*1e-18)
% 
% disp('cellular osmolarity (initial condition)')
% disp(sum(x_c(4:9,1)) + P.chi_C./x_c(3,1)*1e18)
% 
% disp('cellular osmolarity (steady state)')
% disp(sum(y_c(4:9,:)) + P.chi_C./y_c(3,:)*1e18)
% 
% disp('Lumenal osmolarity (initial condition)')
% disp(sum(x_l(1:6,:)) + P.phi_A) 
% disp('Lumenal osmolarity (steady state)')
% disp(sum(y_l(1:6,:)) + P.phi_A) 
% disp('interstitium osmolarity')
% disp(sum(Int) + P.phi_B)

[CellPos, I, IntPos] = form_position_vector(seg_prop);
% CellPos = zeros(1,length(cell_prop));
% for i = 1:length(cell_prop)
%     CellPos(i) = cell_prop{i}.centroid(3);
% end
% [CellPos,I] = sort(CellPos);
% IntPos = lumen_prop.segment(1:end-1);

figure
subplot(2,5,1)
plot(CellPos, y_c(1,I),'--')
hold on
plot(CellPos, y_c(2,I),'--')
hold off
legend('V_A','V_B')
ylabel('mV')
title('Membrane Potential')
subplot(2,5,6)
plot(CellPos, y_c(3,I),'-')
%legend('w_C')
ylabel('um3')
title('Cell Volume')
subplot(2,5,2)
plot(CellPos, y_c(4,I),'-')
hold on
plot(CellPos, y_c(5,I),'-')
hold off
legend('Na_C','K_C')
ylabel('mM')
title('Cellular Concentration')
subplot(2,5,3)
plot(CellPos, y_c(6,I),'-')
hold on
plot(CellPos, y_c(7,I),'-')
hold off
legend('Cl_C','HCO_C')
ylabel('mM')
title('Cellular Concentration')
subplot(2,5,4)
plot(CellPos, -log10(y_c(8,I)*1e-3),'-')
title('Cellular pH')
subplot(2,5,5)
plot(CellPos, y_c(9,I),'-')
legend('CO_C')
ylabel('mM')
ylim([1.27,1.29])
title('Cellular Concentration')

subplot(2,5,7)
plot(IntPos, y_l(1,:))
hold on
plot(IntPos, y_l(2,:))
hold off
legend('Na_A','K_A')
ylabel('mM')
title('Lumenal Concentration')
subplot(2,5,8)
plot(IntPos, y_l(3,:))
hold on
plot(IntPos, y_l(4,:))
hold off
legend('Cl_A','HCO_A')
ylabel('mM')
title('Lumenal Concentration')
subplot(2,5,9)
plot(IntPos, -log10(y_l(5,:)*1e-3))
title('Lumenal pH')
subplot(2,5,10)
plot(IntPos, y_l(6,:))
legend('CO_A')
ylabel('mM')
ylim([1.27,1.29])
title('Lumenal Concentration')

% a(27:40) = [1:0.5:7.5]+a(27:40)
% a(35:40) = [1:0.5:3.5]+a(35:40)
% a(38:40) = [1:0.5:2]+a(38:40)
% b=a/50.03
% figure
% plot(b, y_l(1,:),'LineWidth',2)
% plot(b, y_l(2,:),'LineWidth',2)
% plot(b, y_l(3,:),'LineWidth',2)
% plot(b, y_l(1,:),'LineWidth',2)
% hold on
% plot(b, y_l(2,:),'LineWidth',2)
% plot(b, y_l(3,:),'LineWidth',2)
% legend('Na','K','Cl')
% ylim([0,150])
% ylabel('Concentration of Ions (mM)')
% title('Parotid Gland Duct Model')
% 
% yy_c = y(:,[1 : n_c*9]);
% yy_l = y(:,[n_c*9+1:end]);
% 
% figure
% subplot(1,2,1)
% plot(t,yy_c(:,1),'--')
% hold on
% plot(t,yy_c(:,2),'--')
% plot(t,yy_c(:,3)/10)
% plot(t,yy_c(:,[4,5,6,7]))
% plot(t,-log10(yy_c(:,8)*1e-3))
% plot(t,yy_c(:,9))
% hold off
% legend('V_A','V_B','w_C/10','Na_C','K_C','Cl_C','HCO_C','pH','CO_C')
% title('Cellular Concentration')
% subplot(1,2,2)
% plot(t,yy_l(:,[1,2,3,4]))
% hold on
% plot(t,-log10(yy_l(:,5)*1e-3))
% plot(t,yy_l(:,6))
% hold off
% legend('Na_A','K_A','Cl_A','HCO_A','pH','CO_A')
% title('Lumenal Concentration')