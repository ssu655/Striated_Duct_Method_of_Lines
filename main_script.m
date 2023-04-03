%%
% Author: Shan Su
% Date: May 4th 2021

% This is the main script for running and setting up the parotid gland
% striated duct cell fluid transport model using mesh file as generated
% from experimental scan of mouse parotid gland tissue. The cell geometry
% and lumen geometry are all extracted from the ply mesh files.

% In this script, user can specify the volumetric flow rates of the 
% primary saliva and the discretisation interval length of the lumenal 
% compartment in um. 

% clear

%% Model input setup

L_int = 1; % um length of lumen discretisation interval
PSflow = 7*11.91; % um3/s volumetric primary saliva flow rate

fields = {'Na'; 'K'; 'Cl'; 'HCO'; 'H'; 'CO'};
Int = [140.2; 5.3; 102.6; 24.7; 1000*10^(-7.35); 1.28]; % concentration of interstitium
% PS = [143.5; 5.2; 114.5; 34.2; 1000*10^(-7.35); 1.28];  % concentration of Primary Saliva
PS = [136.95; 6.8; 115.3; 28.47; 7.7260e-05; 1.28];  % concentration of Primary Saliva
CIC = [17; 140; 22; 75; 1000*10^(-7.35); 1.28];  % cellular initial concentration
LIC = [143.5; 5.2; 114.5; 34.2; 1000*10^(-7.35); 1.28]; % lumenal initial concentration
Conc = struct;
Conc.Int = cell2struct(num2cell(Int),fields);
Conc.PS = cell2struct(num2cell(PS),fields);
Conc.CIC = cell2struct(num2cell(CIC),fields);
Conc.LIC = cell2struct(num2cell(LIC),fields);

%% Read mesh file and process raw mesh data

[cell_prop, lumen_prop] = process_mesh_info(L_int);
%% Use full 3D model resolution

s_cell_prop = cell_prop;
s_lumen_prop = lumen_prop;

%% Simplifying the model to one ID and one SD compartment

% [s_cell_prop, s_lumen_prop] = simplify_mesh(cell_prop, lumen_prop);

%% Parameter structure setup

[P_i, P_s] = get_parameters(Conc,PSflow);
s_cell_prop = scale_parameters(s_cell_prop, P_i, P_s);

% Initial condition setup

x = setup_IC(Conc, s_cell_prop, s_lumen_prop);

% f_ODE_noMass(1,x,P,cell_prop,lumen_prop,0);

% Setup and solve the ODE

% f_ODE(1,x,P,cell_prop,lumen_prop,1);
tspan = [0,20000];
P = P_s;

% ===========================================
% % run the version of ODE with mass matrix
% tic
% options = odeset('Mass', @(t,x) mass(t,x,cell_prop,lumen_prop), 'MassSingular', 'yes');
% [t,y] = ode15s(@(t,y) f_ODE(t,y,P,cell_prop,lumen_prop,0), tspan, x, options);
% toc

% ==========================================
% run the version of ODE without mass matrix
tic
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t,y] = ode15s(@(t,y) f_ODE_noMass(t,y,P,s_cell_prop,s_lumen_prop,0,0,0), tspan, x);
toc

% f_ODE_noMass(1,y(end,:),P,cell_prop,lumen_prop,1,0,0);

% Sanity checks

n_c = length(s_cell_prop);
x_c = reshape(x(1 : n_c*9),9,[]); %[9, n_c]
x_l = reshape(x(1+n_c*9 : end),6,[]); %[6, n_l]
y_c = reshape(y(end,1 : n_c*9),9,[]); %[9, n_c]
y_l = reshape(y(end,1+n_c*9 : end),6,[]); %[6, n_l]

% Electroneutrality checks
% 
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

% disp('Electroneutrality check lumen: initial (mM)')
% disp((x_l(2,1)+x_l(3,1)-x_l(4,1)-x_l(5,1)+x_l(6,1)).*s_lumen_prop.disc_volume*1e-18)
% disp('Electroneutrality check lumen: final (mM)')
% disp((y_l(2,:)+y_l(3,:)-y_l(4,:)-y_l(5,:)+y_l(6,:)).*s_lumen_prop.disc_volume*1e-18)

% disp('cellular osmolarity (initial condition)')
% disp(sum(x_c(4:9,1)) + P.chi_C./x_c(3,1)*1e18)

% disp('cellular osmolarity (steady state)')
% disp(sum(y_c(4:9,:)) + P.chi_C./y_c(3,:)*1e18)

% disp('Lumenal osmolarity (initial condition)')
% disp(sum(x_l(1:6,:)) + P.phi_A)
% disp('Lumenal osmolarity (steady state)')
% disp(sum(y_l(1:6,:)) + P.phi_A)
% disp('interstitium osmolarity')
% disp(sum(Int) + P.phi_B)

%% Plotting

IntPos = zeros(1,s_lumen_prop.n_disc);
IntPos(1) = s_lumen_prop.disc_length(1);
for i = 2:s_lumen_prop.n_disc
    out = s_lumen_prop.disc_out_Vec(i);
    IntPos(i) = s_lumen_prop.disc_length(i) + IntPos(out);
end
max_length = max(IntPos);
IntPos = max_length - IntPos;
% IntPos = IntPos(1:58);
% y_l = y_l(:,1:58);

CellPos = zeros(1,length(s_cell_prop)); % [100, 50]
CellType = zeros(2, length(s_cell_prop));
for i = 1:length(s_cell_prop)
    CellPos(i) = s_cell_prop{i}.mean_dist;
    if s_cell_prop{i}.type == "I"
        CellType(:,i) = [1,0];
    else
        CellType(:,i) = [0,1];
    end
end
[CellPos,I] = sort(CellPos); % [50, 100]
CellPos = max_length - CellPos;


figure
subplot(3,2,1)
plot(CellPos, y_c(1,I),'.')
hold on
plot(CellPos, y_c(2,I),'.')
hold off
legend('V_A','V_B')
ylabel('mV')
title('Membrane Potential')
subplot(3,2,2)
w = y_c(3,I);
plot(CellPos(find(CellType(1,I))), w(find(CellType(1,I))),'.')
hold on
plot(CellPos(find(CellType(2,I))), w(find(CellType(2,I))),'.')
hold off
legend('ID', 'SD')
ylabel('\mum^3')
title('Cell Volumn')
subplot(3,2,3)
plot(CellPos, y_c(4,I),'.')
hold on
plot(CellPos, y_c(5,I),'.')
plot(CellPos, y_c(6,I),'.')
plot(CellPos, y_c(7,I),'.')
hold off
legend('Na_C','K_C','Cl_C','HCO_C')
ylabel('mM')
title('Cellular Concentration')
subplot(3,2,4)
w = -log10(y_c(8,I)*1e-3);
plot(CellPos(find(CellType(1,I))), w(find(CellType(1,I))),'.')
hold on
plot(CellPos(find(CellType(2,I))), w(find(CellType(2,I))),'.')
hold off
legend('ID', 'SD')
title('Cellular pH')
subplot(3,2,5)
plot(IntPos, y_l(1,:),'.')
hold on
plot(IntPos, y_l(2,:),'.')
plot(IntPos, y_l(3,:),'.')
plot(IntPos, y_l(4,:),'.')
hold off
legend('Na_A','K_A','Cl_A','HCO_A')
ylabel('mM')
xlabel('Dist along duct (\mum)')
title('Lumenal Concentration')
subplot(3,2,6)
plot(IntPos, -log10(y_l(5,:)*1e-3),'.')
xlabel('Dist along duct (\mum)')
title('Lumenal pH')

%% plot formatting for the simplified cell model
% 
% figure
% CellPos = [2,1];
% IntPos = [2,1];
% subplot(3,2,1)
% plot(CellPos, y_c(1,I),'.','MarkerSize',20)
% hold on
% plot(CellPos, y_c(2,I),'.','MarkerSize',20)
% hold off
% legend('V_A','V_B')
% ylabel('mV')
% set(gca,'XTick','');
% xlim([0,3])
% xlabel('ID cell           SD cell')
% title('Membrane Potential')
% subplot(3,2,2)
% w = y_c(3,I);
% plot(CellPos(find(CellType(1,I))), w(find(CellType(1,I))),'.','MarkerSize',20)
% hold on
% plot(CellPos(find(CellType(2,I))), w(find(CellType(2,I))),'.','MarkerSize',20)
% hold off
% legend('ID', 'SD')
% set(gca,'XTick','');
% xlim([0,3])
% xlabel('ID cell           SD cell')
% ylabel('\mum^3')
% title('Cell Volumn')
% subplot(3,2,3)
% plot(CellPos, y_c(4,I),'.','MarkerSize',20)
% hold on
% plot(CellPos, y_c(5,I),'.','MarkerSize',20)
% plot(CellPos, y_c(6,I),'.','MarkerSize',20)
% plot(CellPos, y_c(7,I),'.','MarkerSize',20)
% hold off
% legend('Na_C','K_C','Cl_C','HCO_C')
% ylabel('mM')
% set(gca,'XTick','');
% xlim([0,3])
% xlabel('ID cell           SD cell')
% title('Cellular Concentration')
% subplot(3,2,4)
% w = -log10(y_c(8,I)*1e-3);
% plot(CellPos(find(CellType(1,I))), w(find(CellType(1,I))),'.','MarkerSize',20)
% hold on
% plot(CellPos(find(CellType(2,I))), w(find(CellType(2,I))),'.','MarkerSize',20)
% hold off
% legend('ID', 'SD')
% set(gca,'XTick','');
% xlim([0,3])
% xlabel('ID cell           SD cell')
% title('Cellular pH')
% subplot(3,2,5)
% plot(IntPos, y_l(1,:),'.','MarkerSize',20)
% hold on
% plot(IntPos, y_l(2,:),'.','MarkerSize',20)
% plot(IntPos, y_l(3,:),'.','MarkerSize',20)
% plot(IntPos, y_l(4,:),'.','MarkerSize',20)
% hold off
% legend('Na_A','K_A','Cl_A','HCO_A')
% ylabel('mM')
% set(gca,'XTick','');
% xlim([0,3])
% xlabel('Intercalated Duct    Striated Duct')
% title('Lumenal Concentration')
% subplot(3,2,6)
% plot(IntPos, -log10(y_l(5,:)*1e-3),'.','MarkerSize',20)
% set(gca,'XTick','');
% xlim([0,3])
% xlabel('Intercalated Duct    Striated Duct')
% title('Lumenal pH')
