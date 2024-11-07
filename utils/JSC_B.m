function A = JSC_B(Pa,Rsa)

A = Pa.J_hB + JSA_B(Pa,Rsa);
end


function JSAB = JSA_B(Pa,Rsa)

Qsso3 = VecToso3(Pa.RH*Rsa*Pa.s);
JSAB = (Pa.RH*Rsa)*Pa.J_saQ*(Pa.RH*Rsa).'-Pa.m_sa*(Pa.dso3^2+Pa.dso3*Qsso3+Qsso3*Pa.dso3);

end