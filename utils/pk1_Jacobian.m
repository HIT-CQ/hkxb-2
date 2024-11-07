function Jaco = pk1_Jacobian(x,Pa,th,p,Tau_sc,Tau_sa)

f = x(1:3);
dth = x(4);

pR = p(1:3);
pth = p(4);
h = Pa.dt;
Rsa = [cos(th) 0 sin(th);0 1 0;-sin(th) 0 cos(th)];
J = JSC_B(Pa,Rsa);
fso3 = VecToso3(f);
gam1 = 0.5*(Pa.RH*Rsa*Pa.J_saQ-Pa.m_sa*Pa.dso3*Pa.RH*Rsa*Pa.sso3)*Pa.e2;
gam2 = 0.5*(Pa.RH*Rsa*Pa.e2so3*Pa.J_saQ-Pa.m_sa*Pa.dso3*Pa.RH*Rsa*Pa.e2so3*Pa.sso3)*Pa.e2;
Gam = gam1*dth;
gam1so3 = VecToso3(gam1);
gam2so3 = VecToso3(gam2);
Gamso3 = VecToso3(Gam);
B = (Pa.RH*Rsa)*(Pa.e2so3*Pa.J_saQ-Pa.J_saQ*Pa.e2so3)*(Pa.RH*Rsa).'...
    -Pa.m_sa*( Pa.dso3*VecToso3(Pa.RH*Rsa*Pa.e2so3*Pa.s)+VecToso3(Pa.RH*Rsa*Pa.e2so3*Pa.s)*Pa.dso3 );
JJ = 0.5*trace(B)*eye(3) - B;
Lmd = JJ+dth*gam2so3-gam1so3;
aph1 = h*pR+0.5*(h^2)*Tau_sc;
aph1_bar = aph1-2*Gam;
aph2 = h*pth+0.5*(h^2)*Tau_sa-Pa.J_saQ2*dth;



dE11 = J+0.5*(fso3*J-VecToso3(J*f))-Gamso3+0.5*Gamso3*fso3-fso3*Gamso3-f*Gam.'-(Gam.'*f)*eye(3)-0.5*aph1_bar*f.';
dE12 = fso3*gam1+0.5*fso3^2*gam1-gam1.'*f*f+0.5*(f.'*f)*gam1+2*gam1;
dE21 = -(so3ToVec(Lmd-Lmd.')).'-0.5*(so3ToVec(fso3*Lmd+Lmd.'*fso3)).'-0.5*(so3ToVec(Lmd*fso3+fso3*Lmd.')).'-0.5*aph2*f.';
dE22 = Pa.J_saQ2+0.25*(f.'*f)*Pa.J_saQ2+trace(fso3*gam2so3)+0.5*trace(fso3*fso3*gam2so3);

Jaco = [dE11 dE12;dE21 dE22];
end