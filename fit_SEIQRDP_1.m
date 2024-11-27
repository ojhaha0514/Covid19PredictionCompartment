function [alpha1,beta1,gamma1,delta1,lambda1,kappa1] = fit_SEIQRDP_1(Q,R,D,Npop,E0,I0,time,guess,varargin)
%% Inputparseer
p = inputParser();
p.CaseSensitive = false;
p.addOptional('tolX',1e-5);  %  option for optimset
p.addOptional('tolFun',1e-5);  %  option for optimset
p.addOptional('Display','iter'); % Display option for optimset
p.addOptional('dt',0.1); % time step for the fitting
p.parse(varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%
tolX = p.Results.tolX ;
tolFun = p.Results.tolFun ;
Display  = p.Results.Display ;
dt  = p.Results.dt ;

%% Options for lsqcurvfit
options=optimset('TolX',tolX,'TolFun',tolFun,...
    'MaxFunEvals',1200,'Display',Display);
%% Initial conditions and basic checks

% Write the target input into a matrix
Q(Q<0)=0; % negative values are not possible
R(R<0)=0; % negative values are not possible
D(D<0)=0; % negative values are not possible

if isempty(R)
    warning(' No data available for "Recovered" ')
    input = [Q;D]; % In this aprticular case, Q is actually the number of active + recovered cases
else
    intermediate = circshift(Q+R+D, 1);
    intermediate(1) = 0;
    intermediate2 = Q+R+D - intermediate;
    input = intermediate2(2:end);
    %input = [Q;R;D];
    %input = Q+R+D;
end

if size(time,1)>size(time,2) && size(time,2)==1,    time = time';end
if size(time,1)>1 && size(time,2)>1,  error('Time should be a vector');end

%% Definition of the new, refined, time vector for the numerical solution
fs = 1./dt;
tTarget = round(datenum(time-time(1))*fs)/fs % Number of days with one decimal
t = tTarget(1):dt:tTarget(end); % oversample to ensure that the algorithm converges

%% Main fitting

modelFun1 = @SEIQRDP_for_fitting; % transform a nested function into anonymous function

ub = [1 1]; % upper bound of the parameters
lb = [0 0]; % lower bound of the parameters
% call Lsqcurvefit
[Coeff] = lsqcurvefit(@(para,t) modelFun1(para,t),...
    guess([2,4]),tTarget(:)',input,lb,ub,options);


%% Write the fitted coeff in the outputs
beta1 = abs(Coeff(1));

alpha1 = 1/14;
gamma1 = 1/2.1;
delta1 = abs(Coeff(2));
lambda1 = 1/25;
kappa1 = lambda1/99;

%% nested functions

    function [output] = SEIQRDP_for_fitting(para,t0)
        
        % I simply rename the inputs
        beta = abs(para(1));
        
        alpha = 1/14;
        gamma = 1/2.1;
        delta = abs(para(2));
        lambda = 1/25;
        kappa = lambda/99;
        
        %% Initial conditions
        N = numel(t);
        Y = zeros(7,N); %  There are seven different states
        Y(2,1) = E0;
        Y(3,1) = I0;
        Y(4,1) = Q(1);
        if ~isempty(R)
            Y(5,1) = R(1);
            Y(1,1) = Npop-Q(1)-R(1)-D(1)-E0-I0;
        else
            Y(1,1) = Npop-Q(1)-D(1)-E0-I0;
        end
        Y(6,1) = D(1);
        
        if round(sum(Y(:,1))-Npop)~=0
            error(['the sum must be zero because the total population',...
                ' (including the deads) is assumed constant']);
        end
        %%
        modelFun = @(Y,A,F) A*Y + F;
        
        % Very large recovery rate should not occur but can lead to
        % numerical errors.
        if lambda>10, warning('lambda is abnormally high'); end
        
        % ODE resolution
        for ii=1:N-1
            A = getA(alpha,gamma,delta,lambda,kappa);
            SI = Y(1,ii)*Y(3,ii);
            F = zeros(7,1);
            F(1:2,1) = [-beta/Npop;beta/Npop].*SI;
            Y(:,ii+1) = RK4(modelFun,Y(:,ii),A,F,dt);
        end
        
        Q1 = Y(4,1:N);
        R1 = Y(5,1:N);
        D1 = Y(6,1:N);
        
        Q1 = interp1(t,Q1,t0);
        R1 = interp1(t,R1,t0);
        D1 = interp1(t,D1,t0);
        if ~isempty(R)
            intermediate3 = circshift(Q1+R1+D1, 1);
            intermediate3(1) = 0;
            intermediate4 = Q1+R1+D1 - intermediate3;
            output = intermediate4(2:end);
            %output = ([Q1;R1;D1]);
            %output = Q1+R1+D1;
            
        else
            output = ([Q1+R1;D1]);
        end
        
    end
    function [A] = getA(alpha,gamma,delta,lambda,kappa)
        %  [A] = getA(alpha,gamma,delta,lambda,kappa) computes the matrix A
        %  that is found in: dY/dt = A*Y + F
        %
        %   Inputs:
        %   alpha: scalar [1x1]: protection rate
        %   beta: scalar [1x1]: infection rate
        %   gamma: scalar [1x1]: Inverse of the average latent time
        %   delta: scalar [1x1]: rate of people entering in quarantine
        %   lambda: scalar [1x1]: cure rate
        %   kappa: scalar [1x1]: mortality rate
        %   Output:
        %   A: matrix: [7x7]
        
        A = zeros(7);
        % S
        A(1,1) = -alpha;
        % E
        A(2,2) = -gamma;
        % I
        A(3,2:3) = [gamma,-delta];
        % Q
        A(4,3:4) = [delta,-kappa-lambda];
        % R
        A(5,4) = lambda;
        % D
        A(6,4) = kappa;
        % P
        A(7,1) = alpha;
        
    end
    function [Y] = RK4(Fun,Y,A,F,dt)
        % NUmerical trick: the parameters are assumed constant between
        % two time steps.
        
        % Runge-Kutta of order 4
        k_1 = Fun(Y,A,F);
        k_2 = Fun(Y+0.5*dt*k_1,A,F);
        k_3 = Fun(Y+0.5*dt*k_2,A,F);
        k_4 = Fun(Y+k_3*dt,A,F);
        % output
        Y = Y + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dt;
    end
end

