function [Q,R,D,T] = getMultipleWaves_2(guess,Npop,time,Confirmed,Recovered,Deaths,Vaccination,tStart1,tIndex,pIndex,fit,est_day,varargin)
%% varargin
Confirmed = Confirmed';
Deaths = Deaths';

%% Inputparseer
p = inputParser();
p.CaseSensitive = false;
p.addOptional('Q0',0.1*Confirmed(1));
p.addOptional('E0',0.2*Confirmed(1));
p.addOptional('I0',0.1*Confirmed(1));
p.parse(varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%
Q0 = p.Results.Q0 ; % initial number of active cases
E0 = p.Results.E0 ; % Initial number of exposed cases. Unknown but unlikely to be zero.
I0 = p.Results.I0 ; % Initial number of infectious cases. Unknown but unlikely to be zero.





%% Remove unecessary data
Confirmed(time<tStart1) = [];
Recovered(time<tStart1) = [];
Deaths(time<tStart1) = [];
time(time<tStart1) = [];

% %  Time for first wave
% indT = [];
% 
% indT(1) = find(time>=tStart1 & time<tIndex(1));
% 
% for jj = 2:length(tIndex)
%     indT(jj) = find(time>=tIndex(jj-1) & time<tIndex(jj));
% end
% 
% indT(length(tIndex)+1) = find(time>=tIndex(end) & time<=tEnd);



%% Simulate first wave
% Initial conditions
R0 = 0.9*Confirmed(find(time==tIndex(1)));
D0 = Deaths(find(time==tIndex(1)));
V0 = Vaccination(find(time==tIndex(1)))*0.875;
U0 = Vaccination(find(time==tIndex(1)))*0.125;
P0 = 0;

E = [];
I = [];
Q = [];
R = [];
D = [];
V = [];
U = [];
P = [];
T = [];


for ii=1:length(tIndex)-1
    indT = find(time>=tIndex(ii) & time<=tIndex(ii+1));
    time1 = time(indT);
    indP = find(pIndex~=0);
    if ~isempty(indP)
        start_delta = tIndex(indP(1));
    else
        start_delta = 0;
    end
    p = pIndex(ii);
    
    [E1,I1,Q1,R1,D1,V1,U1,P1,T1,para] = computeWave(Confirmed(indT),Recovered(indT),...
    Deaths(indT),Vaccination(indT),E0,I0,Q0,R0,D0,V0,U0,P0,time1,time1(1),time1(end),guess,p,start_delta,fit);
  
    parameter = para;

    E = [E, E1];
    I = [I, I1];
    Q = [Q, Q1];
    R = [R, R1];
    D = [D, D1];
    V = [V, V1];
    U = [U, U1];
    P = [P, P1];
    T = [T, T1];
    
    E0 = E(end);
    I0 = I(end);
    Q0 = Q(end);
    R0 = R(end);
    D0 = D(end);
    V0 = V(end);
    U0 = U(end);
    P0 = P(end);

end

time_pred = time(end):time(end)+est_day;
dt0 =1/24;
newT_pred = time(end):dt0:time(end)+est_day;
Vaccination_pred1 = repelem(median(Vaccination(end-6:end)), est_day+1);
Vaccination_pred2 = repelem(median(Vaccination(end-6:end))*1.5, est_day+1);
Vaccination_pred3 = repelem(median(Vaccination(end-6:end))*0.5, est_day+1);
N0 = numel(newT_pred);
t_pred = [0:N0-1].*dt0;

[~,E2,I2,Q2,R2,D2,V2,P2] = SEIQRDVUP_2(parameter(1),parameter(2),parameter(3),parameter(4),parameter(5),...
            Vaccination_pred1,Npop,E0,I0,Q0,R0,D0,V0,U0,P0,time_pred,t_pred,p,start_delta,fit);
[~,E3,I3,Q3,R3,D3,V3,P3] = SEIQRDVUP_2(0.75*parameter(1),parameter(2),parameter(3),parameter(4),parameter(5),...
            Vaccination_pred1,Npop,E0,I0,Q0,R0,D0,V0,U0,P0,time_pred,t_pred,p,start_delta,fit);
[~,E4,I4,Q4,R4,D4,V4,P4] = SEIQRDVUP_2(1.25*parameter(1),parameter(2),parameter(3),parameter(4),parameter(5),...
            Vaccination_pred1,Npop,E0,I0,Q0,R0,D0,V0,U0,P0,time_pred,t_pred,p,start_delta,fit);
        
        
        E = [E, E2, E3, E4];
        I = [I, I2, I3, I4];
        Q = [Q, Q2, Q3, Q4];
        R = [R, R2, R3, R4];
        D = [D, D2, D3, D4];
        V = [V, V2, V3, V4];
        P = [P, P2, P3, P4];
        T = [T, newT_pred, newT_pred, newT_pred];

%% Nested functions

    function [E,I,Q,R,D,V,U,P,newT,para] = computeWave(Confirmed,Recovered,Deaths,Vaccination,...
            E0,I0,Q0,R0,D0,V0,U0,P0,time,tStart,tEnd,guess,p,start_delta,fit)
        
        % Parameter estimation with the lsqcurvefit function
        [mu1,kappa1,alpha1,gamma1,omega1] = ...
            fit_SEIQRDVUP_2(Confirmed,Recovered,Deaths,Vaccination,Npop,...
            E0,I0,V0,U0,P0,time,guess,p,start_delta,fit,'Display','off');
        
        para = [mu1, kappa1, alpha1, gamma1, omega1];
        
        dt = 1/24; % time step
        newT = tStart:dt:tEnd;
        N = numel(newT);
        t = [0:N-1].*dt;
        [~,E,I,Q,R,D,V,U,P] = SEIQRDVUP_2(mu1,kappa1,alpha1,gamma1,omega1,...
            Vaccination,Npop,E0,I0,Q0,R0,D0,V0,U0,P0,time,t,p,start_delta,fit);
    end

end

