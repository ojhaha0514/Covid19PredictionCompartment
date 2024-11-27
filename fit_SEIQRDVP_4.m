function [mu1,kappa1,alpha1,gamma1,omega1] = fit_SEIQRDVP_4(C,Vaccination,Npop,N_age,E0,I0,Q0,R0,D0,V0,P0,time,guess,pIndex,start_delta,fit,beta,varargin)
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
%dt  = p.Results.dt ;
dt = 1/24;
%% Options for lsqcurvfit
opt = optimset('Display','off');
%% Initial conditions and basic checks

% Write the target input into a matrix
C(C<0)=0; % negative values are not possible
input = sum(C,2)


%% Definition of the new, refined, time vector for the numerical solution
% fs = 1./dt;
% tTarget = round(datenum(time-time(1))*fs)/fs; % Number of days with one decimal
% t = tTarget(1):dt:tTarget(end) % oversample to ensure that the algorithm converges

t = 0:dt:(length(time)-1);
Vacc = repelem(Vaccination/24, 24);

%% Main fitting
% fit_proportion = variants_proportion();
modelFun1 = @SEIQRDVP_for_fitting; % transform a nested function into anonymous function

ub = 10; % upper bound of the parameters
lb = 0; % lower bound of the parameters
% call Lsqcurvefit
[Coeff] = lsqcurvefit(@(para,t) modelFun1(para,t),...
    guess(1),(1:length(C(:,1))-1),input,lb,ub,opt);


%% Write the fitted coeff in the outputs
mu1 = abs(Coeff);
alpha1 = 1/6;

omega1 = 1/43;
kappa1 = 1/2.1;
gamma1 = 1/25;


%% nested functions

    function [output] = SEIQRDVP_for_fitting(para,t0)
        
        
        mu = abs(para);
        alpha = 1/6;
        
        omega = 1/43;
        kappa = 1/2.1;
        gamma = 1/25;
        f = 0.01;
        
        %% Initial conditions
        N = numel(t);
        Y = zeros(8,N); %  There are seven different states
    
        Y(1,1) = N_age(1)-C(1,1)-V0-(E0+P0)*N_age(1)/sum(N_age);
        Y(2,1) = E0*N_age(1)/sum(N_age);
        Y(3,1) = I0*N_age(1)/sum(N_age);
        Y(4,1) = Q0*N_age(1)/sum(N_age);
        Y(5,1) = R0*N_age(1)/sum(N_age);
        Y(6,1) = D0*N_age(1)/sum(N_age);
        Y(7,1) = V0; 
        Y(8,1) = P0*N_age(1)/sum(N_age);
        
        for ii = 2:9 
           Y(1,1,ii) = N_age(ii)-C(1,ii)-V0-(E0+P0)*N_age(ii)/sum(N_age);
           Y(2,1,ii) = E0*N_age(ii)/sum(N_age);
           Y(3,1,ii) = I0*N_age(ii)/sum(N_age);
           Y(4,1,ii) = Q0*N_age(ii)/sum(N_age);
           Y(5,1,ii) = R0*N_age(ii)/sum(N_age);
           Y(6,1,ii) = D0*N_age(ii)/sum(N_age);
           Y(7,1,ii) = V0;
           Y(8,1,ii) = P0*N_age(ii)/sum(N_age);
        end
            
        %%
        x0 = datenum(time(1))-datenum(datetime(start_delta));
        
        beta = ones(9,9);
        
        for ii=1:100
            for jj=1:9
                s0 = Y(1,ii,jj);
                e0 = Y(2,ii,jj);
                i0 = beta(1,jj)*Y(3,ii,1)+beta(2,jj)*Y(3,ii,2)+beta(3,jj)*Y(3,ii,3)+...
                    beta(4,jj)*Y(3,ii,4)+beta(5,jj)*Y(3,ii,5)+beta(6,jj)*Y(3,ii,6)+...
                    beta(7,jj)*Y(3,ii,7)+beta(8,jj)*Y(3,ii,8)+beta(9,jj)*Y(3,ii,9);
                q0 = Y(4,ii,jj);
                r0 = Y(5,ii,jj);
                d0 = Y(6,ii,jj);
                v0 = Y(7,ii,jj);
                p0 = Y(8,ii,jj);
                
                if isnan(i0)
                   break 
                end
                
                if pIndex
                    
                    fun = @(x,xdata)exp(x(1)*xdata)./(x(2)+exp(x(1)*xdata));
                    
                    p_var = fun(fit,(x0+(ii-1)/24)/100);
                else
                    p_var = 0;
                end
                
                e_var = 0.1*(1-p_var)+(0.1+0.05)*p_var;
                
                s1 = -(1+p_var*0.97)*(mu)*s0*i0/Npop-Vacc(ii);
                e1 = (1+p_var*0.97)*(mu)*(s0+e_var*v0)*i0/Npop - kappa*e0;
                i1 = kappa*e0 - alpha*i0;
                q1 = alpha*i0 - gamma*q0;
                r1 = (1-f)*gamma*q0;
                d1 = f*gamma*q0;
                v1 = Vacc(ii) - (1+p_var*0.97)*(mu)*e_var*v0*i0/Npop - omega*v0;
                p1 = omega*v0;
                
                if pIndex
                    p_var = fun(fit,(x0+(ii-1)/24+dt/2)/100);
                else
                    p_var = 0;
                end
                
                e_var = 0.1*(1-p_var)+(0.1+0.05)*p_var;
                
                s2 = -(1+p_var*0.97)*(mu)*(s0+dt*s1/2)*(i0+dt*i1/2)/Npop-Vacc(ii);
                e2 = (1+p_var*0.97)*(mu)*(s0+e_var*v0+dt*(s1+e_var*v1)/2)*(i0+dt*i1/2)/Npop - kappa*(e0+dt*e1/2);
                i2 = kappa*(e0+dt*e1/2) - alpha*(i0+dt*i1/2);
                q2 = alpha*(i0+dt*i1/2) - gamma*(q0+dt*q1/2);
                r2 = (1-f)*gamma*(q0+dt*q1/2);
                d2 = f*gamma*(q0+dt*q1/2);
                v2 = Vacc(ii) - (1+p_var*0.97)*(mu)*e_var*(v0+dt*v1/2)*(i0+dt*i1/2)/Npop - omega*(v0+dt*v1/2);
                p2 = omega*(v0+dt*v1/2);
                
                s3 = -(1+p_var*0.97)*(mu)*(s0+dt*s2/2)*(i0+dt*i2/2)/Npop-Vacc(ii);
                e3 = (1+p_var*0.97)*(mu)*(s0+e_var*v0+dt*(s2+e_var*v2)/2)*(i0+dt*i2/2)/Npop - kappa*(e0+dt*e2/2);
                i3 = kappa*(e0+dt*e2/2) - alpha*(i0+dt*i2/2);
                q3 = alpha*(i0+dt*i2/2) - gamma*(q0+dt*q2/2);
                r3 = (1-f)*gamma*(q0+dt*q2/2);
                d3 = f*gamma*(q0+dt*q2/2);
                v3 = Vacc(ii) - (1+p_var*0.97)*(mu)*e_var*(v0+dt*v2/2)*(i0+dt*i2/2)/Npop - omega*(v0+dt*v2/2);
                p3 = omega*(v0+dt*v2/2);
                
                if pIndex
                    p_var = fun(fit,(x0+(ii-1)/24+dt)/100);
                else
                    p_var = 0;
                end
                
                e_var = 0.1*(1-p_var)+(0.1+0.05)*p_var;
                
                s4 = -(1+p_var*0.97)*(mu)*(s0+dt*s3)*(i0+dt*i3)/Npop-Vacc(ii);
                e4 = (1+p_var*0.97)*(mu)*(s0+e_var*v0+dt*(s3+e_var*v3))*(i0+dt*i3)/Npop - kappa*(e0+dt*e3);
                i4 = kappa*(e0+dt*e3) - alpha*(i0+dt*i3);
                q4 = alpha*(i0+dt*i3) - gamma*(q0+dt*q3);
                r4 = (1-f)*gamma*(q0+dt*q3);
                d4 = f*gamma*(q0+dt*q3);
                v4 = Vacc(ii) - (1+p_var*0.97)*(mu)*e_var*(v0+dt*v3)*(i0+dt*i3)/Npop - omega*(v0+dt*v3);
                p4 = omega*(v0+dt*v3);
                
                Y(1,ii+1,jj) = s0+dt*(s1+2*s2+2*s3+s4)/6;
                Y(2,ii+1,jj) = e0+dt*(e1+2*e2+2*e3+e4)/6;
                Y(3,ii+1,jj) = i0+dt*(i1+2*i2+2*i3+i4)/6;
                Y(4,ii+1,jj) = q0+dt*(q1+2*q2+2*q3+q4)/6;
                Y(5,ii+1,jj) = r0+dt*(r1+2*r2+2*r3+r4)/6;
                Y(6,ii+1,jj) = d0+dt*(d1+2*d2+2*d3+d4)/6;
                Y(7,ii+1,jj) = v0+dt*(v1+2*v2+2*v3+v4)/6;
                Y(8,ii+1,jj) = p0+dt*(p1+2*p2+2*p3+p4)/6;
                
            end
        end
        Y
%         Q1 = Y(4,1:N);
%         R1 = Y(5,1:N);
%         D1 = Y(6,1:N);
%         
%         Q1 = interp1(t,Q1,t0);
%         R1 = interp1(t,R1,t0);
%         D1 = interp1(t,D1,t0);

        output = [];
        for jj = 1:9
            Q1 = Y(4,1:24:N,jj);
            R1 = Y(5,1:24:N,jj);
            D1 = Y(6,1:24:N,jj);
            output = [output (Q1+R1+D1)'];
        end
        output = sum(output,2)
    end
end

