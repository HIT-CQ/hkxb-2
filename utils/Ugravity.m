function [FUsc,TauUsc,TauUsa] = Ugravity(Pa,r,R,th)

rnorm3 = (norm(r))^3;
rnorm5 = (norm(r))^5;
muA3 = Pa.mu/rnorm3;
muA5 = Pa.mu/rnorm5;
Rsa = [cos(th) 0 sin(th);0 1 0;-sin(th) 0 cos(th)];
a = (R*Pa.RH*Rsa).'*r+(Pa.RH*Rsa).'*Pa.d;

FUsc = -muA3*Pa.m_sc*(r);
TauUsc = 3*muA5*cross(R.'*r,JSC_B(Pa,Rsa)*R.'*r);
TauUsa = Pa.e2.'*3*muA5*cross(a,Pa.J_saQ*a);

end

