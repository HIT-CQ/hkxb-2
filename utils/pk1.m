function E = pk1(x,Pa,th,p,Tau_sc,Tau_sa)

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
B = (Pa.RH*Rsa)*(Pa.e2so3*Pa.J_saQ-Pa.J_saQ*Pa.e2so3)*(Pa.RH*Rsa).'...
    -Pa.m_sa*( Pa.dso3*VecToso3(Pa.RH*Rsa*Pa.e2so3*Pa.s)+VecToso3(Pa.RH*Rsa*Pa.e2so3*Pa.s)*Pa.dso3 );
JJ = 0.5*trace(B)*eye(3) - B;
Lmd = JJ+dth*VecToso3(gam2)-VecToso3(gam1);
aph1 = h*pR  +0.5*(h^2)*Tau_sc;
aph1_bar = aph1-2*Gam;
aph2 = h*pth +0.5*(h^2)*Tau_sa-Pa.J_saQ2*dth;



E1 = J*f+0.5*fso3*(J*f)+fso3*Gam+0.5*(fso3^2)*Gam...
     - Gam.'*f*f - 0.25*(f.'*f)*aph1_bar - aph1_bar;
E2 = trace(fso3*Lmd)+0.5*trace(fso3*fso3*Lmd)-0.25*(f.'*f)*aph2 - aph2;

E = [E1;E2];
end