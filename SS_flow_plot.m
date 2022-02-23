
flow = [1:5:90];
lumen_data = zeros(6,length(flow));
for i = 1:length(flow)
    lumen_data(:,i) = dynamic_script(flow(i), 0, P, y,s_cell_prop,s_lumen_prop);
end

%% %
load('MouseParotid.mat')
gland = MouseParotid; % RabbitSM;

% load('MouseSM.mat')
% gland = MouseSM; % RabbitSM;

scale = min(gland.K(:,1));
xmax = 22;
flow = flow*0.25;
col = get(gca,'colororder');

figure
subplot(1,2,1)
plot(flow, lumen_data(1,:),'-','LineWidth',2,'color',col(1,:))
hold on
plot(gland.Na(:,1)/scale,gland.Na(:,2),"o",'color',col(1,:))
plot(flow, lumen_data(2,:),'-','LineWidth',2,'color',col(2,:))
plot(gland.K(:,1)/scale,gland.K(:,2),"o",'color',col(2,:))
hold off
xlim([0,xmax])
legend("Simulation - Na", "Experiment - Na","Simulation - K", "Experiment - K")
% legend("Experiment - Na", "Experiment - K")
% xlabel("Flow Rate (\mu l/min g)")
xlabel("Multiples of unstimulated saliva flow rate")
ylabel("Concentration (mM)")
subplot(1,2,2)
plot(flow, lumen_data(3,:),'-','LineWidth',2,'color',col(1,:))
hold on
plot(gland.Cl(:,1)/scale,gland.Cl(:,2),"o",'color',col(1,:))
plot(flow, lumen_data(4,:),'-','LineWidth',2,'color',col(2,:))
plot(gland.HCO(:,1)/scale,gland.HCO(:,2),"o",'color',col(2,:))
hold off
xlim([0,xmax])
xlabel("Multiples of unstimulated saliva flow rate")
legend("Simulation - Cl", "Experiment - Cl","Simulation - HCO", "Experiment - HCO")
% legend("Experiment - Cl", "Experiment - HCO")
