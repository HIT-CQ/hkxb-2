% Time
Pa.dt = 0.1;
Pa.end_time = 1200;
Pa.t = 0:Pa.dt:Pa.end_time;

% Stacked Plate Satellite
Pa.m_h  = 175.220;  % kg
Pa.m_sa = 12.06;
Pa.m_sc = Pa.m_h+Pa.m_sa;
Pa.J_h  = [14.07 -1.60 -0.79;-1.60 71.36 -0.32;-0.79 -0.32 84.23];     % kg，m2 
Pa.J_sa = [2.04 -0.036 -0.01;-0.036 5.27 -0.00049;-0.01 -0.00049 7.31];
Pa.RH = [-1 0 0;0 cos(pi/4) sin(pi/4);0 -sin(pi/4) cos(pi/4)];
Pa.c_B = [0.00986;0.00965;0.13662]; % m
Pa.d = [-1.185;0.577;0.203];
Pa.s = [1;0;0];         % m
Pa.Ax = 3.1;
Pa.Ay = 0.55;
Pa.Az = 5;
Pa.CoP = [0.54; 0.20; -0.082];
Pa.CD = 2.2;

% Common
Pa.e1 = [1;0;0];
Pa.e2 = [0;1;0];
Pa.e3 = [0;0;1];

% Intermediate variable
Pa.c_Bso3 = VecToso3(Pa.c_B);
Pa.sso3 = VecToso3(Pa.s);
Pa.dso3 = VecToso3(Pa.d);
Pa.e2so3 = VecToso3(Pa.e2);
Pa.J_hB  = Pa.J_h - Pa.m_h*(Pa.c_Bso3)^2;
Pa.J_saQ = Pa.J_sa - Pa.m_sa*(Pa.sso3)^2;
Pa.J_saQ2 = Pa.e2.'*Pa.J_saQ*Pa.e2;
Pa.K = 300; % (N，m)/rad
Pa.C = 0;  % (N，m，s)/rad

% Ini value
Ini.f_guess = [0;0;0];
Ini.delta_th_guess = 0;
Ini.th0 = 0;
Ini.R0 = eye(3);
Ini.Rsa0 = [cos(Ini.th0) 0 sin(Ini.th0);0 1 0;-sin(Ini.th0) 0 cos(Ini.th0)];
Ini.w0 = [0;0;0];
Ini.dth0 = 0;
Ini.r0 = [6948137;0;0];
Ini.v0 = [0;4969.1;5716.3];
Ini.p0 = mass(Pa,Ini.Rsa0)*[Ini.w0;Ini.dth0];

% Earth
Pa.mu = 3.986e14;
Pa.r_earth = 6378.137e3;
Pa.we = [0;0.7292e-4;0];
Pa.rho = 1.742e-13; % 570km
Pa.J2 = 0.00108263;