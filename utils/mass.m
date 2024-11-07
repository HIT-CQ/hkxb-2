  function MM = mass(Pa,Rsa)

M11 = JSC_B(Pa,Rsa);
M12 = (Pa.RH*Rsa*Pa.J_saQ-Pa.m_sa*Pa.dso3*Pa.RH*Rsa*Pa.sso3)*Pa.e2;

M21 = M12.';
M22 = Pa.J_saQ2;

MM = [M11  M12;
      M21  M22];
end