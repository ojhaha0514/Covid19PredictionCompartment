function [S,E,I,Q,R,D,V,P] = SEIQRDVP_3(mu,kappa,alpha,gamma,omega,Vaccination,Npop,E0,I0,Q0,R0,D0,V0,P0,time,t,pIndex,start_delta,fit)
%% Initial conditions
mu;
N = numel(t);
Vacc = repelem(Vaccination/24, 24);
dt = median(diff(t));
Y = zeros(8,N); %  There are seven different states
        
Y(1,1) = Npop-Q0-R0-D0-E0-I0-V0-P0;
Y(2,1) = E0;
Y(3,1) = I0;
Y(4,1) = Q0;
Y(5,1) = R0;
Y(6,1) = D0;
Y(7,1) = V0;
Y(8,1) = P0;

if round(sum(Y(:,1))-Npop)~=0
    error(['the sum must be zero because the total population',...
        ' (including the deads) is assumed constant']);
end
%% Computes the nine states

% ODE resolution

f = 0.01;
% x0 = datenum(time(1))-datenum(start_delta);

for ii=1:N-1
    s0 = Y(1,ii);
    e0 = Y(2,ii);
    i0 = Y(3,ii);
    q0 = Y(4,ii);
    r0 = Y(5,ii);
    d0 = Y(6,ii);
    v0 = Y(7,ii);
    p0 = Y(8,ii);
    
    if pIndex
        
        fun = @(x,xdata)exp(x(1)*xdata)./(x(2)+exp(x(1)*xdata));
        p_var = fun(fit,(x0+(ii-1)/24)/100);
    else
        p_var = 0;
    end
    
    e_var = 0.22*(1-p_var)+(0.1+0.05)*p_var;
    
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
    
    e_var = 0.22*(1-p_var)+(0.1+0.05)*p_var;
    
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
    
    e_var = 0.22*(1-p_var)+(0.1+0.05)*p_var;
    
    s4 = -(1+p_var*0.97)*(mu)*(s0+dt*s3)*(i0+dt*i3)/Npop-Vacc(ii);
    e4 = (1+p_var*0.97)*(mu)*(s0+e_var*v0+dt*(s3+e_var*v3))*(i0+dt*i3)/Npop - kappa*(e0+dt*e3);
    i4 = kappa*(e0+dt*e3) - alpha*(i0+dt*i3);
    q4 = alpha*(i0+dt*i3) - gamma*(q0+dt*q3);
    r4 = (1-f)*gamma*(q0+dt*q3);
    d4 = f*gamma*(q0+dt*q3);
    v4 = Vacc(ii) - (1+p_var*0.97)*(mu)*e_var*(v0+dt*v3)*(i0+dt*i3)/Npop - omega*(v0+dt*v3);
    p4 = omega*(v0+dt*v3);
    
    Y(1,ii+1) = s0+dt*(s1+2*s2+2*s3+s4)/6;
    Y(2,ii+1) = e0+dt*(e1+2*e2+2*e3+e4)/6;
    Y(3,ii+1) = i0+dt*(i1+2*i2+2*i3+i4)/6;
    Y(4,ii+1) = q0+dt*(q1+2*q2+2*q3+q4)/6;
    Y(5,ii+1) = r0+dt*(r1+2*r2+2*r3+r4)/6;
    Y(6,ii+1) = d0+dt*(d1+2*d2+2*d3+d4)/6;
    Y(7,ii+1) = v0+dt*(v1+2*v2+2*v3+v4)/6;
    Y(8,ii+1) = p0+dt*(p1+2*p2+2*p3+p4)/6;
    
end

% Y = round(Y);
%% Write the outputs
S = Y(1,1:N);
E = Y(2,1:N);
I = Y(3,1:N);
Q = Y(4,1:N);
R = Y(5,1:N);
D = Y(6,1:N);
V = Y(7,1:N);
P = Y(8,1:N);


%% Nested functions
   
end


