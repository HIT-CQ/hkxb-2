
%{

 Nonlinear Lie Group variational computational model  
of plate-type broadband satellites in LEO aerodynamic flow       
                                                                        
Author: Qian Cao                                                       
Date: Oct. 28th, 2024    
                                               
%}

clc
close all
clear

%% Add function folder path
addpath('utils');

%% Load Parameter
Parameter;

%% Solver

% ode45 Solution
disp('----------------------------------------');
disp('----------------------------------------');
disp('                ODE Method:             ');
disp('----------------------------------------');
disp('----------------------------------------');

x0 = [Ini.r0;reshape(Ini.R0,9,1);Ini.th0;Ini.v0;Ini.p0];
[t,x] = ode45(@(t,x) odeDyn(t,x,Pa),[0:Pa.dt:Pa.end_time],x0);
[~,whist1,dthhist1] = cellfun(@(t,x) odeDyn(t,x.',Pa), num2cell(t), num2cell(x,2),'uni',0);
% [t, x] = RK4(@(t, x) odeDyn(t, x, Pa), [0 Pa.end_time], x0, Pa.dt);
% [~,whist1,dthhist1] = cellfun(@(t,x) odeDyn(t,x.',Pa), num2cell(t), num2cell(x,2),'uni',0);
rhist1 = x(1:end,1:3).';
Rhist1_column = x(1:end,4:12).';
for j = 1:length(t)
    Rhist1(:,3*j-2:3*j) = reshape(Rhist1_column(:,j), 3, 3);
    error_Rhist1(:,j) = norm(Rhist1(:,3*j-2:3*j).'*Rhist1(:,3*j-2:3*j)-eye(3));
end
thhist1 = x(1:end,13).';
thhist1 = rad2deg(thhist1);
whist1 = reshape(cell2mat(whist1),3,length(t));
whist1 = rad2deg(whist1);
dthhist1 = reshape(cell2mat(dthhist1),1,length(t));
dthhist1 = rad2deg(dthhist1);


% My Solution
disp('------------------------------------------');
disp('------------------------------------------');
disp('Lie Group Variational Computational Model:');
disp('------------------------------------------');
disp('------------------------------------------');

[t,rhist,Fhist,whist,dthhist,Rhist,error_Rhist,thhist,phist] = satLVI(Pa,Ini);
whist = rad2deg(whist);
thhist = rad2deg(thhist);
dthhist = rad2deg(dthhist);

%% Test Plot
figure('name','Orbit')
plot(t, rhist(1,:), 'k-', t, rhist(2,:), 'b--',t, rhist(3,:))
legend('r_1','r_2','r_3')
xlabel('Time/s')
ylabel('Orbit Position/m')
title('Time-Position Curve')

figure('name','Vibration-Angle')
plot(t, thhist, 'k-')
legend('$\theta$','Interpreter', 'latex')
xlabel('Time/s')
ylabel('Vibration-Angle $\theta$/rad','Interpreter', 'latex')
title('Time-Vibration Angle Curve')

figure('name','ODE-Vibration-Angle')
plot(t, thhist1, 'k-')
legend('$\theta$','Interpreter', 'latex')
xlabel('Time/s')
ylabel('Vibration-Angle $\theta$/rad','Interpreter', 'latex')
title('ODE-Time-Vibration Angle Curve')

figure('name','Vibration-Angle Velocity')
plot(t(1:end-1), dthhist, 'k-')
legend('$\dot{\theta}$','Interpreter', 'latex')
xlabel('Time/s')
ylabel('Vibration-Angle Velocity $\theta$ (rad/s)','Interpreter', 'latex')
title('Time-Vibration Angle Velocity Curve')

figure('name','ODE-Vibration-Angle Velocity')
plot(t, dthhist1, 'k-')
legend('$\dot{\theta}$','Interpreter', 'latex')
xlabel('Time/s')
ylabel('Vibration-Angle Velocity $\theta$ (rad/s)','Interpreter', 'latex')
title('ODE-Time-Vibration Angle Velocity Curve')

figure('name','Angular Velocity')
plot(t(1:end-1), whist(1,:), 'k-', t(1:end-1), whist(2,:), 'b--',t(1:end-1), whist(3,:),'k:')
legend('w_1','w_2','w_3')
xlabel('Time/s')
ylabel('Angular Velocity (deg/s)')
title('Time-Position Curve')

figure('name','ODE-Angular Velocity')
plot(t, whist1(1,:), 'k-', t, whist1(2,:), 'b--',t, whist1(3,:),'k:')
legend('w_1','w_2','w_3')
xlabel('Time/s')
ylabel('Angular Velocity (deg/s)')
title('ODE-Time-Position Curve')

figure('name','Orthogonality Error')
plot(t, error_Rhist, 'k-')
legend('$||R^{T} R - I_3 ||$','Interpreter', 'latex')
xlabel('Time/s')
ylabel('Orthogonality Error $||R^{T} R - I_3||$','Interpreter', 'latex')
title('Orthogonality Error Curve')

figure('name','ODE-Orthogonality Error')
plot(t, error_Rhist1, 'k-')
legend('$||R^{T} R - I_3 ||$','Interpreter', 'latex')
xlabel('Time/s')
ylabel('Orthogonality Error $||R^{T} R - I_3||$','Interpreter', 'latex')
title('ODE-Orthogonality Error Curve')

%% Remove function folder path
rmpath('utils');