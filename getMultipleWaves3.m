function [Q,R,D,T] = getMultipleWaves3(guess,Npop,time,Confirmed,Recovered,Deaths,tStart1,tIndex,varargin)
%% varargin

Active = Confirmed-Recovered-Deaths;
Active(Active<0) = 0; % No negative number possible

%% Inputparseer
p = inputParser();
p.CaseSensitive = false;
p.addOptional('Q0',Active(1));
p.addOptional('E0',0.3*Active(1));
p.addOptional('I0',5*Active(1));
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
Active = Confirmed-Recovered-Deaths;
Active(Active<0) = 0; % No negative number possible

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
R0 = Recovered(find(time==tIndex(1)));
D0 = Deaths(find(time==tIndex(1)));

E = [];
I = [];
Q = [];
R = [];
D = [];
T = [];


for ii=1:length(tIndex)-1
    indT = find(time>=tIndex(ii) & time<=tIndex(ii+1));
    time1 = time(indT);
    [E1,I1,Q1,R1,D1,T1] = computeWave(Active(indT),Recovered(indT),...
    Deaths(indT),E0,I0,Q0,R0,D0,time1,time1(1),time1(end),guess);
  


    E = [E, E1];
    I = [I, I1];
    Q = [Q, Q1];
    R = [R, R1];
    D = [D, D1];
    T = [T, T1];
    
    E0 = E(end);
    I0 = I(end);
    Q0 = Q(end);
    R0 = R(end);
    D0 = D(end);

end

%% Nested functions

    function [E,I,Q,R,D,newT] = computeWave(Active,Recovered,Deaths,E0,I0,Q0,R0,D0,time,tStart,tEnd,guess)
        
        % Parameter estimation with the lsqcurvefit function
        [alpha1,beta1,gamma1,delta1,Lambda1,Kappa1,lambdaFun,kappaFun] = ...
            fit_SEIQRDP(Active,Recovered,Deaths,Npop,E0,I0,time,guess,'Display','off');
        
        dt = 1/24; % time step
        newT = tStart:dt:tEnd;
        N = numel(newT);
        t = [0:N-1].*dt;
        [~,E,I,Q,R,D,~] = SEIQRDP(alpha1,beta1,gamma1,delta1,Lambda1,Kappa1,...
            Npop,E0,I0,Q0,R0,D0,t,lambdaFun,kappaFun);
    end

    function [rmse] = RMSE(y1,y2)
        
        rmse = sqrt(nanmean((y1(:)-y2(:)).^2));
        
    end

end

