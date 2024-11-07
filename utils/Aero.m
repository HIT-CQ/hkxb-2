function tau_aer = Aero(r,v,R,Pa)

vref = v - VecToso3(Pa.we)*r;
vref_I_unit = vref/norm(vref);
Areas = [Pa.Ax;Pa.Ay;Pa.Az;-Pa.Ax;-Pa.Ay;-Pa.Az];
l = vref_I_unit.'*R.'*[eye(3) -eye(3)];
l_bar = max(l,0);
F_aer = -0.5*Pa.rho*Pa.CD*(vref.'*vref)*l_bar*Areas;
tau_aer = F_aer*(R.'*VecToso3(vref_I_unit)*Pa.CoP);

end

