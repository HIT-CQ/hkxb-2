function F = Cay2F(f)

fso3 = VecToso3(f);
F = eye(3) + 4/(4+(f.'*f))*(fso3 + 0.5*(fso3^2));

end