% load("result_high.mat")
function lumen_data = dynamic_script(flow, displ, P, y,cell_prop,lumen_prop)
% time_series.time = [0:0.1:500]; % 14.43*7*ones(1,5001);% 
% time_series.Q =  14.43*7*(flow - 1)./(1+exp(-0.1.*(time_series.time-200)))+14.43*7;
% time_series.Na = P.ConP.Na*ones(1,5001);
% time_series.Cl = P.ConP.Cl*ones(1,5001);
% time_series.K = P.ConP.K*ones(1,5001);

load("result_high.mat")
plot(time_series.time,time_series.Q)

step = 0.5;
tspan = [0:step:1000];%[0:step:400,401:25000];
% tspan = [0,1000];
x = y(end,:);

tic
[t,z] = ode15s(@(t,z) f_ODE_noMass(t,z,P,cell_prop,lumen_prop,0,1,time_series), tspan, x);
toc

%% plotting dynamic result
if displ
% load("low_stim.mat")
cell_no = 20;
cell_no = cell_no - 1;
n_c = length(cell_prop);
loc_disc = find(cell_prop{cell_no}.api_area_discs~=0);
yy_c = z(:,[cell_no*9+1 : cell_no*9+9]);
yy_l = z(:,[n_c*9+loc_disc(1)*6+1:n_c*9+loc_disc(1)*6+6]);

figure
subplot(3,2,1)
plot(t, yy_c(:,1),'LineWidth',1)
hold on
plot(t, yy_c(:,2),'LineWidth',1)
hold off
legend('V_A','V_B')
ylabel('mV')
ylim([-80,0])
title('Membrane Potentials')
subplot(3,2,2)
plot(t, yy_c(:,3),'LineWidth',1)
ylabel('\mu m^3')
string1 = strcat('Cell Volume, cell at :',num2str((cell_no+1)*10), '\mu m');
title(string1)
subplot(3,2,3)
plot(t, yy_c(:,4),'LineWidth',1)
hold on
plot(t, yy_c(:,5),'LineWidth',1)
plot(t, yy_c(:,6),'LineWidth',1)
plot(t, yy_c(:,7),'LineWidth',1)
hold off
legend('Na_C','K_C','Cl_C','HCO_C')
ylabel('mM')
title('Cellular Concentrations')
subplot(3,2,4)
plot(t, -log10(yy_c(:,8)*1e-3),'-','LineWidth',1)
title('Cellular pH')
subplot(3,2,5)
plot(t, yy_l(:,1),'LineWidth',1)
hold on
plot(t, yy_l(:,2),'LineWidth',1)
plot(t, yy_l(:,3),'LineWidth',1)
plot(t, yy_l(:,4),'LineWidth',1)
hold off
legend('Na_A','K_A','Cl_A','HCO_A')
ylabel('mM')
xlabel('time (s)')
title('Luminal Concentrations')
subplot(3,2,6)
plot(t, -log10(yy_l(:,5)*1e-3),'LineWidth',1)
xlabel('time (s)')
title('Luminal pH')

% plot whole duct at a fixed time point
time = 400;% second
time_ind = time/step +1;
yyy_c = reshape(z(time_ind,1 : n_c*9),9,[]); %[9, n_c]
yyy_l = reshape(z(time_ind,1+n_c*9 : end),6,[]); %[6, n_l]

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
for i = 1:length(cell_prop)
    CellPos(i) = cell_prop{i}.mean_dist;
end
[CellPos,I] = sort(CellPos);
CellPos = max_length - CellPos;

figure
subplot(3,2,1)
plot(CellPos, yyy_c([1,2],I),'LineWidth',1)
legend('V_A','V_B')
ylabel('mV')
title('Membrane Potential')
subplot(3,2,2)
plot(CellPos, yyy_c(3,I),'-','LineWidth',1)
ylabel('\mu m^3')
title('Cell Volume')
subplot(3,2,3)
plot(CellPos, yyy_c([4,5,6,7],I),'LineWidth',1)
legend('Na_C','K_C','Cl_C','HCO_C')
ylabel('mM')
title('Cellular Concentration')
subplot(3,2,4)
plot(CellPos, -log10(yyy_c(8,I)*1e-3),'LineWidth',1)
title('Cellular pH')
subplot(3,2,5)
plot(IntPos, yyy_l([1,2,3,4],:),'LineWidth',1)
legend('Na_A','K_A','Cl_A','HCO_A')
ylabel('mM')
xlabel('Duct Length (\mum)')
title('Luminal Concentration')
subplot(3,2,6)
plot(IntPos, -log10(yyy_l(5,:)*1e-3),'LineWidth',1)
xlabel('Duct Length (\mum)')
title('Luminal pH')

% figure
% subplot(3,1,1)
% plot(t, yy_c(:,3),'LineWidth',1)
% ylabel('\mu m^3')
% string1 = strcat('Cell Volume, cell at :',num2str((cell_no+1)*10), '\mu m');
% title(string1)
% subplot(3,1,2)
% plot(t, yy_l(:,1),'LineWidth',1)
% hold on
% plot(t, yy_l(:,2),'LineWidth',1)
% plot(t, yy_l(:,3),'LineWidth',1)
% plot(t, yy_l(:,4),'LineWidth',1)
% hold off
% legend('Na_A','K_A','Cl_A','HCO_A')
% ylabel('mM')
% xlabel('time (s)')
% title('Luminal Concentrations')
% subplot(3,1,3)
% plot(t, -log10(yy_l(:,5)*1e-3),'LineWidth',1)
% xlabel('time (s)')
% title('Luminal pH')

% record whole duct info at a fixed time point
time = 500;% second
time_ind = time/step +1;
yyy_c = reshape(z(time_ind,1 : n_c*9),9,[]); %[9, n_c]
yyy_l = reshape(z(time_ind,1+n_c*9 : end),6,[]); %[6, n_l]

end

lumen_data = yyy_l(:,end);
lumen_data = z;
end

