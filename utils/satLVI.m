function [t,rhist,Fhist,whist,dthhist,Rhist,error_Rhist,thhist,phist] = satLVI(Pa,Ini)

% Initialize Variables
h = Pa.dt;
t = Pa.t;
rhist = zeros(3, length(t)-1);
Fhist = zeros(3, 3*(length(t)-1));
vhist = zeros(3, length(t)-1);
whist = zeros(3, length(t)-1);
dthhist = zeros(1, length(t)-1);
Rhist = zeros(3, 3*length(t));
error_Rhist = zeros(1,length(t));
thhist = zeros(1, length(t));
phist = zeros(4, length(t));
prhist = zeros(3,length(t));
Force_U = zeros(3,length(t));
TauUsc = zeros(3,length(t));
TauUsa = zeros(1,length(t));
Tau_sc = zeros(3,length(t));
tau_aer = zeros(3,length(t));
Tau_c = zeros(3,length(t));
Tau_sa = zeros(1,length(t));

rhist(:,1) = Ini.r0;
Rhist(:,1:3) = eye(3);
error_Rhist(:,1) = 0;
thhist(:,1) = 0;
vhist(:,1) = Ini.v0;
dthhist(:,1) = 0;
phist(:,1) = mass(Pa,Ini.Rsa0)*[Ini.w0;Ini.dth0];
prhist(:,1) = Pa.m_sc*Ini.v0;

% Print Indicator
next = 10;
lastPrint = 0;

for k = 1:length(t)-1

    % Print Simulation Progress
    if t(k) == lastPrint
        disp(['Time = ',num2str(t(k)),' out of ',num2str(t(end))])
        lastPrint = lastPrint + next;
    end
    
    % Force, Torque for satellite
    [Force_U(:,k),TauUsc(:,k),TauUsa(:,k)] = Ugravity(Pa,rhist(:,k),Rhist(:,(3*k-2):3*k),thhist(:,k));
    Tau_c(:,k) = Tau_Control(t(k));
    tau_aer(:,k) = Aero(rhist(:,k),vhist(:,k),Rhist(:,(3*k-2):3*k),Pa);
    Tau_sc(:,k) = TauUsc(:,k) + Tau_c(:,k) + tau_aer(:,k);
    Tau_sa(:,k) = TauUsa(:,k) - Pa.K*thhist(:,k) - Pa.C*dthhist(:,k);
    
    % calculate (delta_r)
    delta_r = (prhist(:,k)+h/2*Force_U(:,k))*h/Pa.m_sc;
    
    % Use Newton's method to calculate (f, delta_th)
    fun = @(x) pk1(x,Pa,thhist(:,k),phist(:,k),Tau_sc(:,k),Tau_sa(:,k));
    Jacobi = @(x) pk1_Jacobian(x,Pa,thhist(:,k),phist(:,k),Tau_sc(:,k),Tau_sa(:,k));
    x0 = [Ini.f_guess;Ini.delta_th_guess];
    x = rootn_newton(fun,Jacobi,x0);
    f = x(1:3);
    delta_th = x(4);
    Fhist(:,(3*k-2):3*k) = Cay2F(f);
    whist(:,k) = 1/h*so3ToVec(Fhist(:,(3*k-2):3*k)-eye(3));
    vhist(:,k) = delta_r/h;
    dthhist(:,k) = delta_th/h;

    % Update k+1 States
    rhist(:,k+1) = rhist(:,k) + delta_r;
    Rhist(:,3*k+1:3*(k+1)) = Rhist(:,(3*k-2):3*k)*Fhist(:,(3*k-2):3*k);
    thhist(:,k+1) = thhist(:,k) + delta_th;
    error_Rhist(:,k+1) = norm(Rhist(:,3*k+1:3*(k+1)).'*Rhist(:,3*k+1:3*(k+1))-eye(3));
    
    % Update k+1 generalized momentum
    [Force_U(:,k+1),TauUsc(:,k+1),TauUsa(:,k+1)] = Ugravity(Pa,rhist(:,k+1),Rhist(:,3*k+1:3*(k+1)),thhist(:,k+1));
    Tau_c(:,k+1) = Tau_Control(t(k+1));
    tau_aer(:,k+1) = Aero(rhist(:,k+1),vhist(:,k),Rhist(:,(3*k+1):3*k+3),Pa);
    Tau_sc(:,k+1) = TauUsc(:,k+1) + Tau_c(:,k+1) + tau_aer(:,k+1);
    Tau_sa(:,k+1) = TauUsa(:,k+1) - Pa.K*thhist(:,k+1) - Pa.C*dthhist(:,k);
    
    prhist(:,k+1) = prhist(:,k) + h/2*(Force_U(:,k)+Force_U(:,k+1));
    phist(:,k+1) = pk2(x,Pa,thhist(:,k),Tau_sc(:,k+1),Tau_sa(:,k+1));   
end

end