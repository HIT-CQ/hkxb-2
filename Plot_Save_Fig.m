%% LineStyle  1.Straight Line  2.Dashed Line 3.Dot Line  4.Dotted Line
Lstyle = {'-','--',':','-.'};

%% Color
Cstyle ={ [0 0 0],  [0.5 0.5 0.5],     [1,0,0],     ...       % 1.Black 2. Gray 3. Red
          [0 0 1],  [0 0.4471 0.7412], [0.6 0.8 1], ...       % 4.Deep Blue 5. Sky Blue 6. Light Blue
          [0 1 1],  [0 1 0],           [0.494 0.184 0.556]};  % 7.Cyan 8, Green 9, Purple
      
%% Coordinate Axis
xaxis = '时间 (s)';
yaxis = {'ODE45正交误差','LGVI正交误差',...
         '$\omega_x$ (deg/s)','$\omega_y$ (deg/s)','$\omega_z$ (deg/s)',...
         '$\theta$ (deg)','$\dot{\theta}$ (deg/s)'};

%% Legend
mylegend = {'LGVI','ODE45'}; 

%% Another time
t2 = t(1:end-1);

%% Plot
% Orthogonality Error-ODE
figure('name','Orthogonality Error_ODE')
ha = tightPlots(1, 1, 14, [1 0.618], [0,0], [1.3 0.6], [1.2 0.5], 'centimeters');

axes(ha(1));
plot(t,error_Rhist1(1,:),'color',Cstyle{1},'Linewidth',1,'linestyle',Lstyle{1});
xlabel(xaxis,'FontSize',12,'Color','k');
% ylabel(yaxis{1},'FontSize',12,'Color','k','fontname', 'Monospaced','Interpreter','latex');
ylabel(yaxis{1},'FontSize',12,'Color','k','fontname', 'Monospaced');
% set(legend('FontSize',12,'Box', 'off','fontname', 'Monospaced','orientation','horizontal','location',[0.65 0.8 0.01 0.05]),'Interpreter','latex');

print(gcf, 'Fig1.png', '-dpng', '-r600');
% print(gcf, 'Fig1.tif', '-dtiff', '-r600', '-opengl');
% print(gcf, 'Fig1.pdf', '-dpdf', '-painters', '-loose');
% print(gcf, 'Fig1.eps', '-depsc2', '-painters', '-loose');

% Orthogonality Error-LGVI
figure('name','Orthogonality Error_LGVI')
ha = tightPlots(1, 1, 14, [1 0.618], [0,0], [1.3 0.6], [1.3 0.5], 'centimeters');

axes(ha(1));
plot(t,error_Rhist(1,:),'color',Cstyle{1},'Linewidth',1,'linestyle',Lstyle{1});
xlabel(xaxis,'FontSize',12,'Color','k');
ylabel(yaxis{2},'FontSize',12,'Color','k','fontname', 'Monospaced');

print(gcf, 'Fig2.png', '-dpng', '-r600');
% print(gcf, 'Fig2.tif', '-dtiff', '-r600', '-opengl');
% print(gcf, 'Fig2.pdf', '-dpdf', '-painters', '-loose');
% print(gcf, 'Fig2.eps', '-depsc2', '-painters', '-loose');

% Angular Velocity-x
figure('name','wx')
ha = tightPlots(1, 1, 14, [1 0.618], [0,0], [1.3 0.6], [1.35 0.5], 'centimeters');

axes(ha(1));
plot(t2,whist(1,:),'color',Cstyle{1},'Linewidth',1,'linestyle',Lstyle{1},'displayname',mylegend{1});
hold on
plot(t,whist1(1,:),'color',Cstyle{3},'Linewidth',1,'linestyle',Lstyle{1},'displayname',mylegend{2});
xlabel(xaxis,'FontSize',12,'Color','k');
ylabel(yaxis{3},'FontSize',12,'Color','k','fontname', 'Monospaced','Interpreter','latex');
set(legend('orientation','horizontal','location',[0.7 0.2 0.01 0.05]), 'FontSize',12,'Box', 'off');

axes('Position',[0.3 0.35 0.5 0.5]);
plot(t2,whist(1,:),'color',Cstyle{1},'Linewidth',1,'linestyle',Lstyle{1});
hold on
plot(t,whist1(1,:),'color',Cstyle{3},'Linewidth',1,'linestyle',Lstyle{1});
axis([40 200 74 78]);

print(gcf, 'Fig3.png', '-dpng', '-r600');
% print(gcf, 'Fig3.tif', '-dtiff', '-r600', '-opengl');
% print(gcf, 'Fig3.pdf', '-dpdf', '-painters', '-loose');
% print(gcf, 'Fig3.eps', '-depsc2', '-painters', '-loose');

% Angular Velocity-y
figure('name','wy')
ha = tightPlots(1, 1, 14, [1 0.618], [0,0], [1.3 0.6], [1.35 0.5], 'centimeters');

axes(ha(1));
plot(t2,whist(2,:),'color',Cstyle{1},'Linewidth',1,'linestyle',Lstyle{1},'displayname',mylegend{1});
hold on
plot(t,whist1(2,:),'color',Cstyle{3},'Linewidth',1,'linestyle',Lstyle{1},'displayname',mylegend{2});
xlabel(xaxis,'FontSize',12,'Color','k');
ylabel(yaxis{4},'FontSize',12,'Color','k','fontname', 'Monospaced','Interpreter','latex');
set(legend('orientation','horizontal','location',[0.7 0.27 0.01 0.05]), 'FontSize',12,'Box', 'off');

axes('Position',[0.3 0.4 0.5 0.5]);
plot(t2,whist(2,:),'color',Cstyle{1},'Linewidth',1,'linestyle',Lstyle{1});
hold on
plot(t,whist1(2,:),'color',Cstyle{3},'Linewidth',1,'linestyle',Lstyle{1});
axis([40 200 -10 -8.5]);

print(gcf, 'Fig4.png', '-dpng', '-r600');
% print(gcf, 'Fig4.tif', '-dtiff', '-r600', '-opengl');
% print(gcf, 'Fig4.pdf', '-dpdf', '-painters', '-loose');
% print(gcf, 'Fig4.eps', '-depsc2', '-painters', '-loose');

% Angular Velocity-z
figure('name','wz')
ha = tightPlots(1, 1, 14, [1 0.618], [0,0], [1.3 0.6], [1.35 0.5], 'centimeters');

axes(ha(1));
plot(t2,whist(3,:),'color',Cstyle{1},'Linewidth',1,'linestyle',Lstyle{1},'displayname',mylegend{1});
hold on
plot(t,whist1(3,:),'color',Cstyle{3},'Linewidth',1,'linestyle',Lstyle{1},'displayname',mylegend{2});
xlabel(xaxis,'FontSize',12,'Color','k');
ylabel(yaxis{5},'FontSize',12,'Color','k','fontname', 'Monospaced','Interpreter','latex');
set(legend('orientation','vertical','location',[0.85 0.82 0.01 0.05]), 'FontSize',12,'Box', 'off');

axes('Position',[0.3 0.5 0.4 0.4]);
plot(t2,whist(3,:),'color',Cstyle{1},'Linewidth',1,'linestyle',Lstyle{1});
hold on
plot(t,whist1(3,:),'color',Cstyle{3},'Linewidth',1,'linestyle',Lstyle{1});
axis([40 200 -4 -2]);

print(gcf, 'Fig5.png', '-dpng', '-r600');
% print(gcf, 'Fig5.tif', '-dtiff', '-r600', '-opengl');
% print(gcf, 'Fig5.pdf', '-dpdf', '-painters', '-loose');
% print(gcf, 'Fig5.eps', '-depsc2', '-painters', '-loose');

% Vibration Angle-ODE45
figure('name','Vibration Angle-ODE45')
ha = tightPlots(1, 1, 14, [1 0.618], [0,0], [1.3 0.6], [1.35 0.5], 'centimeters');

axes(ha(1));
plot(t,thhist1(1,:),'color',Cstyle{1},'Linewidth',1,'linestyle',Lstyle{1});
xlabel(xaxis,'FontSize',12,'Color','k');
ylabel(yaxis{6},'FontSize',12,'Color','k','fontname', 'Monospaced','Interpreter','latex');

print(gcf, 'Fig6.png', '-dpng', '-r600');
% print(gcf, 'Fig6.tif', '-dtiff', '-r600', '-opengl');
% print(gcf, 'Fig6.pdf', '-dpdf', '-painters', '-loose');
% print(gcf, 'Fig6.eps', '-depsc2', '-painters', '-loose');


% Vibration Angle-LGVI
figure('name','Vibration Angle-LGVI')
ha = tightPlots(1, 1, 14, [1 0.618], [0,0], [1.3 0.6], [1.35 0.5], 'centimeters');

axes(ha(1));
plot(t,thhist(1,:),'color',Cstyle{1},'Linewidth',1,'linestyle',Lstyle{1});
xlabel(xaxis,'FontSize',12,'Color','k');
ylabel(yaxis{6},'FontSize',12,'Color','k','fontname', 'Monospaced','Interpreter','latex');

print(gcf, 'Fig7.png', '-dpng', '-r600');
% print(gcf, 'Fig7.tif', '-dtiff', '-r600', '-opengl');
% print(gcf, 'Fig7.pdf', '-dpdf', '-painters', '-loose');
% print(gcf, 'Fig7.eps', '-depsc2', '-painters', '-loose');