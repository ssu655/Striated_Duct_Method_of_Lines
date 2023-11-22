% load("result_high.mat")
% function lumen_data = dynamic_script(flow, displ, P, y,cell_prop,lumen_prop)
% flow = 3;
% time_series.time = [0:0.1:500]; % 14.43*7*ones(1,5001);% 
% time_series.Q = 11.91*7*(flow - 1)./(1+exp(-0.1.*(time_series.time-100)))+11.91*7;
% time_series.Na = P.ConP.Na*ones(1,5001); %14.43*7
% time_series.Cl = P.ConP.Cl*ones(1,5001);
% time_series.HCO = P.ConP.HCO*ones(1,5001);
% time_series.H = P.ConP.H*ones(1,5001);
% time_series.K = P.ConP.K*ones(1,5001);

% load("Acinus PDE Results\result_bicarb_smooth_VPLC0.002.mat")
load("result_bicarb_VPLC0.004.mat")
time_series.Q = time_series.Q*7;

% for i = 7800:length(time_series.time)
%     if i+14<length(time_series.time)
%         window = i-14:i+14;
%     else
%         window = i-14:length(time_series.time);
%     end
%     time_series.Q(i) = mean(time_series.Q(window));
%     time_series.Na(i) = mean(time_series.Na(window));
%     time_series.K(i) = mean(time_series.K(window));
%     time_series.Cl(i) = mean(time_series.Cl(window));
%     time_series.HCO(i) = mean(time_series.HCO(window));
%     time_series.H(i) = mean(time_series.H(window));
% end
% 
% figure (3)
% subplot(3,2,1)
% plot(time_series.time,time_series.Q)
% subplot(3,2,2)
% plot(time_series.time,time_series.Na)
% subplot(3,2,3)
% plot(time_series.time,time_series.K)
% subplot(3,2,4)
% plot(time_series.time,time_series.Cl)
% subplot(3,2,5)
% plot(time_series.time,time_series.HCO)
% subplot(3,2,6)
% plot(time_series.time,time_series.H)

%%
step = 0.1;
tspan1 = [0:0.1:400];%[0:step:400,401:25000];
x = y(end,:);

tic
[t,z1] = ode15s(@(t,z) f_ODE_noMass(t,z,P,s_cell_prop,s_lumen_prop,0,1,time_series), tspan1, x);
toc

tspan2 = [400.1:0.1:800];
x = z1(end,:);

tic
[t,z2] = ode15s(@(t,z) f_ODE_noMass(t,z,P,s_cell_prop,s_lumen_prop,0,1,time_series), tspan2, x);
toc

t = [tspan1,tspan2];
z = [z1;z2];
%% Plot the fluxes at a fixed time point (steady-state with stimulation)

ind = find(t==380); % at time = 400 second

f_ODE_noMass(t(ind), z(ind,:), P,s_cell_prop, s_lumen_prop, 1, 1, time_series);

%% plot whole duct at a fixed time point (steady-state with stimulation)
time = 380;% second
time_ind = time/step +1;
n_c = length(cell_prop);
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
CellType = zeros(2, length(cell_prop));
for i = 1:length(cell_prop)
    CellPos(i) = cell_prop{i}.mean_dist;
    if cell_prop{i}.type == "I"
        CellType(:,i) = [1,0];
    else
        CellType(:,i) = [0,1];
    end
end
[CellPos,I] = sort(CellPos);
CellPos = max_length - CellPos;

figure
ax(1) = subplot(3,2,1);
plot(CellPos, yyy_c([1,2],I),'.','MarkerSize',10)
legend('V_A','V_B','Orientation','horizontal','Location','southoutside')
ylabel('mV')
title('Membrane Potential')
ax(2) = subplot(3,2,2);
w = yyy_c(3,I);
plot(CellPos(find(CellType(1,:))), w(find(CellType(1,:))),'.','MarkerSize',10)
hold on
plot(CellPos(find(CellType(2,:))), w(find(CellType(2,:))),'.','MarkerSize',10)
hold off
legend('ID', 'SD','Orientation','horizontal','Location','southoutside')
ylabel('\mu m^3')
title('Cell Volume')
ax(3) = subplot(3,2,3);
plot(CellPos, yyy_c([4,5,6,7],I),'.','MarkerSize',10)
legend('Na_C','K_C','Cl_C','HCO_C','Orientation','horizontal','Location','southoutside')
ylabel('mM')
title('Cellular Concentration')
ax(4) = subplot(3,2,4);
plot(CellPos, -log10(yyy_c(8,I)*1e-3),'.','MarkerSize',10)
title('Cellular pH')
hLegend = legend('','Location','southoutside');
set(hLegend,'visible','off')
ax(5) = subplot(3,2,5);
plot(IntPos, yyy_l([1,2,3,4],:),'.','MarkerSize',10)
legend('Na_A','K_A','Cl_A','HCO_A','Orientation','horizontal','Location','southoutside')
ylabel('mM')
xlabel('ID entry     SD entry                                SD exit')
title('Local Duct Concentration')
ax(6) = subplot(3,2,6);
plot(IntPos, -log10(yyy_l(5,:)*1e-3),'.','MarkerSize',10)
xlabel('ID entry     SD entry                                SD exit')
title('Local Duct pH')
hLegend = legend('','Location','southoutside');
set(hLegend,'visible','off')

sgtitle('Steady-state duct solution') 
set(gcf,'position',[250,50,800,700])
x_range = [-5,140];

for k = 1:6
    xlim(ax(k),x_range)
    set(ax(k),'xtick',[],'YGrid','on','xlim',x_range)
end
%% plotting single cell dynamic result
% if displ
% load("low_stim.mat")
cell_no = 20;
cell_no = cell_no - 1;
n_c = length(cell_prop);
poc = cell_prop{cell_no}.mean_dist;
loc_disc = find(cell_prop{cell_no}.api_area_discs~=0);
yy_c = z(:,[cell_no*9+1 : cell_no*9+9]);
yy_l = z(:,[n_c*9+loc_disc(1)*6+1:n_c*9+loc_disc(1)*6+6]);

figure
ax(1) = subplot(3,2,1);
plot(t, yy_c(:,1),'LineWidth',2)
hold on
plot(t, yy_c(:,2),'LineWidth',2)
hold off
legend('V_A','V_B','Orientation','horizontal','Location','southoutside')
ylabel('mV')
ylim([-80,0])
title('Membrane Potentials')
ax(2) = subplot(3,2,2);
plot(t, yy_c(:,3),'LineWidth',2)
ylabel('\mum^3')
string1 = strcat('Volume of a cell half way along th SD');
title(string1)
hLegend = legend('','Location','southoutside');
set(hLegend,'visible','off')
ax(3) = subplot(3,2,3);
plot(t, yy_c(:,4),'LineWidth',2)
hold on
plot(t, yy_c(:,5),'LineWidth',2)
plot(t, yy_c(:,6),'LineWidth',2)
plot(t, yy_c(:,7),'LineWidth',2)
hold off
legend('Na_C','K_C','Cl_C','HCO_C','Orientation','horizontal','Location','southoutside')
ylabel('mM')
title('Cellular Concentrations')
ax(4) = subplot(3,2,4);
plot(t, -log10(yy_c(:,8)*1e-3),'-','LineWidth',2)
title('Cellular pH')
hLegend = legend('','Location','southoutside');
set(hLegend,'visible','off')
ax(5) = subplot(3,2,5);
plot(t, yy_l(:,1),'LineWidth',2)
hold on
plot(t, yy_l(:,2),'LineWidth',2)
plot(t, yy_l(:,3),'LineWidth',2)
plot(t, yy_l(:,4),'LineWidth',2)
hold off
legend('Na_A','K_A','Cl_A','HCO_A','Orientation','horizontal','Location','southoutside')
ylabel('mM')
xlabel('time (s)')
title('Local duct Concentrations')
ax(6) = subplot(3,2,6);
plot(t, -log10(yy_l(:,5)*1e-3),'LineWidth',2)
xlabel('time (s)')
hLegend = legend('','Location','southoutside');
set(hLegend,'visible','off')
title('Local duct pH')

sgtitle('Single duct cell') 
set(gcf,'position',[300,150,800,700])
x_range = [0,800];

for k = 1:6
    xlim(ax(k),x_range)
%     set(ax(k),'xtick',[],'YGrid','on','xlim',x_range)
end
%% showing whole duct info across all time
% time = 500;% second
% time_ind = time/step +1;
% yyy_c = reshape(z(time_ind,1 : n_c*9),9,[]); %[9, n_c]
% yyy_l = reshape(z(time_ind,1+n_c*9 : end),6,[]); %[6, n_l]

% t_sample = 1:200:10000;
t_sample = 1:100:5000;
x = 1:2:(length(lumen_prop.disc_length)-1);
pos = zeros(size(x));
tt = t(t_sample);
NAA = zeros(length(x),length(tt));
KA = zeros(length(x),length(tt));
CLA = zeros(length(x),length(tt));

% looping through all the duct discs
for i = 1:length(x)
% cell_no = x(i);
% cell_no = cell_no - 1;
n_c = length(cell_prop);
% yy_c = z(:,[cell_no*9+1 : cell_no*9+9]);
yy_l = z(:,[n_c*9+x(i)*6+1:n_c*9+x(i)*6+6]);
pos(i) = IntPos(x(i));
NAA(i,:) = yy_l(t_sample,1);
KA(i,:) = yy_l(t_sample,2);
CLA(i,:) = yy_l(t_sample,3);
end

figure
subplot(2,2,1)
[ppos,I] = sort(pos);
[X,Y] = meshgrid(ppos,tt);
surf(X,Y,NAA(I,:)')
zlabel('mM')
ylabel('Time (s)')
xlabel('Dist along Duct (\mum)')
title('A : [Na^+]_A')

subplot(2,2,2)
surf(X,Y,KA(I,:)')
zlabel('mM')
ylabel('Time (s)')
xlabel('Dist along Duct (\mum)')
title('B : [K^+]_A')

subplot(2,2,3)
surf(X,Y,CLA(I,:)')
zlabel('mM')
ylabel('Time (s)')
xlabel('Dist along Duct (\mum)')
title('C : [Cl^-]_A')

x = 2:28:length(cell_prop);
x = [90,86,76,58,30,2];
cpos = zeros(size(x));
tt = t(t_sample);
wC = zeros(length(x),length(tt));

for i = 1:length(x)
cell_no = x(i);
cell_no = cell_no - 1;
n_c = length(cell_prop);
cell_prop{x(i)}.type
yy_c = z(:,[cell_no*9+1 : cell_no*9+9]);
cpos(i) = max_length - cell_prop{x(i)}.mean_dist;
wC(i,:) = yy_c(t_sample,3);
end

subplot(2,2,4) 
[ppos,I] = sort(cpos);
[X,Y] = meshgrid(ppos,tt);
surf(X,Y,wC(I,:)')
zlabel('\mum^3')
ylabel('Time (s)')
xlabel('Dist along Duct (\mum)')
xlim([0,150])
title('D : Cell volume')

lumen_data = yyy_l(:,1);
% lumen_data = z;
% end

