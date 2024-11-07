function Tau_C = Tau_Control(t)

if t>=30 && t<=40
    Tau_C = [3;0;0];
else
    Tau_C = [0;0;0];
end

% Tau_C = [0;0;0];

end