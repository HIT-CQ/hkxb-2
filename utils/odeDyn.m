function [dstate,w,dth] = odeDyn(t,state,Pa)
r = state(1:3);
R = reshape(state(4:12), 3, 3);
th = state(13);
v = state(14:16);
pw = state(17:19);
pth = state(20);

%% Intermediate Variables
Rsa = [cos(th) 0 sin(th);0 1 0;-sin(th) 0 cos(th)];
[FUsc,TauUsc,TauUsa] = Ugravity(Pa,r,R,th);
B = (Pa.RH*Rsa)*(Pa.e2so3*Pa.J_saQ-Pa.J_saQ*Pa.e2so3)*(Pa.RH*Rsa).' ...
    - Pa.m_sa*(Pa.dso3*VecToso3(Pa.RH*Rsa*Pa.e2so3*Pa.s)+VecToso3(Pa.RH*Rsa*Pa.e2so3*Pa.s)*Pa.dso3);

dq = mass(Pa,Rsa) \ state(17:20);
w = dq(1:3);


%% Kinematics
dr = v;
dR = R*VecToso3(w);
dth = dq(4);

%% Force & Torque
Tau_sc  = TauUsc + Aero(r,v,R,Pa) + Tau_Control(t);
Tau_sa = TauUsa - Pa.K*th - Pa.C*dth;

%% Dynamics
dv = FUsc/Pa.m_sc;
dpw = - VecToso3(w)*pw + Tau_sc;
dpth = 0.5*w.'*B*w ... 
       + w.'*(Pa.RH*Rsa*Pa.e2so3*Pa.J_saQ-Pa.m_sa*Pa.dso3*Pa.RH*Rsa*Pa.e2so3*Pa.sso3)*Pa.e2*dth ...
       + Tau_sa;

dstate = [dr;reshape(dR,9,1);dth;dv;dpw;dpth];
end