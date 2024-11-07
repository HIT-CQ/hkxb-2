function pNext = pk2(x,Pa,th,Tau_sc,Tau_sa)

f = x(1:3);
dth = x(4);

h = Pa.dt;
Rsa = [cos(th) 0 sin(th);0 1 0;-sin(th) 0 cos(th)];
J = JSC_B(Pa,Rsa);
F = Cay2F(f);
gam1 = 0.5*(Pa.RH*Rsa*Pa.J_saQ-Pa.m_sa*Pa.dso3*Pa.RH*Rsa*Pa.sso3)*Pa.e2;
gam1so3 = VecToso3(gam1);
Gam = gam1*dth;
a = 0.5*trace(J)*eye(3) - J + VecToso3(Gam);
Phi2 = so3ToVec(a*F-F.'*a.');



pR = 1/h*Phi2 + h/2*Tau_sc;
pth = 1/h*Pa.J_saQ2*dth + 1/h*trace((eye(3)-F)*gam1so3) + h/2*Tau_sa;


pNext = [pR;pth];
end