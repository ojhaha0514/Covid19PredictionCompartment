function [S,E,I,Q,R,D,V1,V2,V3,P] = SEIQRDVP_5(mu,kappa,alpha,gamma,omega,Vaccination1,Vaccination2,Vaccination3,Npop,E0,I0,Q0,R0,D0,V10,V20,V30,P0,time,t,pIndex,start_delta,fit)
%% Initial conditions
N = numel(t);
Vacc1 = repelem(Vaccination1/24, 24);
Vacc2 = repelem(Vaccination2/24, 24);
Vacc3 = repelem(Vaccination3/24, 24);
dt = median(diff(t));
Y = zeros(10,N); %  There are seven different states
        
Y(1,1) = Npop-Q0-R0-D0-E0-I0-V10-V20-V30-P0;
Y(2,1) = E0;
Y(3,1) = I0;
Y(4,1) = Q0;
Y(5,1) = R0;
Y(6,1) = D0;
Y(7,1) = V10;
Y(8,1) = V20;
Y(9,1) = V30;
Y(10,1) = P0;

if round(sum(Y(:,1))-Npop)~=0
    error(['the sum must be zero because the total population',...
        ' (including the deads) is assumed constant']);
end
%% Computes the nine states

% ODE resolution

f = 0.01;

x0 = datenum(time(1))-datenum(start_delta);


for ii=1:N-1
    s0 = Y(1,ii);
    e0 = Y(2,ii);
    i0 = Y(3,ii);
    q0 = Y(4,ii);
    r0 = Y(5,ii);
    d0 = Y(6,ii);
    v10 = Y(7,ii);
    v20 = Y(8,ii);
    v30 = Y(9,ii);
    p0 = Y(10,ii);
    
    if pIndex

        fun = @(x,xdata)exp(x(1)*xdata)./(x(2)+exp(x(1)*xdata));

        p_var = fun(fit,(x0+(ii-1)/24)/100);
    else
        p_var = 0;
    end

    e_var1 = 0.25;
    e_var2 = 0.2;
    e_var3 = 0.15;
    
    s1 = -(1+p_var*(pIndex-1))*(mu)*s0*i0/Npop-Vacc1(ii);
    e1 = (1+p_var*(pIndex-1))*(mu)*(s0+(e_var1*v10+e_var2*v20+e_var3*v30))*i0/Npop - kappa*e0;
    i1 = kappa*e0 - alpha*i0;
    q1 = alpha*i0 - gamma*q0;
    r1 = (1-f)*gamma*q0;
    d1 = f*gamma*q0;
    v11 = Vacc1(ii)-Vacc2(ii) - (1+p_var*(pIndex-1))*(mu)*e_var1*v10*i0/Npop;
    v21 = Vacc2(ii)-Vacc3(ii) - (1+p_var*(pIndex-1))*(mu)*e_var2*v20*i0/Npop;
    v31 = Vacc3(ii) - (1+p_var*(pIndex-1))*(mu)*e_var3*v30*i0/Npop - omega*v30;
    p1 = omega*v30;

    %             if pIndex
    %                 p_var = fun(fit,(x0+(ii-1)/24+dt/2)/100);
    %             else
    %                 p_var = 0;
    %             end
    %

    s2 = -(1+p_var*(pIndex-1))*(mu)*(s0+dt*s1/2)*(i0+dt*i1/2)/Npop-Vacc1(ii);
    e2 = (1+p_var*(pIndex-1))*(mu)*(s0+(e_var1*v10+e_var2*v20+e_var3*v30)+...
        dt*(s1+(e_var1*v11+e_var2*v21+e_var3*v31))/2)*(i0+dt*i1/2)/Npop - kappa*(e0+dt*e1/2);
    i2 = kappa*(e0+dt*e1/2) - alpha*(i0+dt*i1/2);
    q2 = alpha*(i0+dt*i1/2) - gamma*(q0+dt*q1/2);
    r2 = (1-f)*gamma*(q0+dt*q1/2);
    d2 = f*gamma*(q0+dt*q1/2);
    v12 = Vacc1(ii) - Vacc2(ii) - (1+p_var*(pIndex-1))*(mu)*e_var1*(v10+dt*v11/2)*(i0+dt*i1/2)/Npop;
    v22 = Vacc2(ii) - Vacc3(ii) - (1+p_var*(pIndex-1))*(mu)*e_var2*(v20+dt*v21/2)*(i0+dt*i1/2)/Npop;
    v32 = Vacc3(ii) - (1+p_var*(pIndex-1))*(mu)*e_var3*(v30+dt*v31/2)*(i0+dt*i1/2)/Npop - omega*(v30+dt*v31/2);
    p2 = omega*(v30+dt*v31/2);

    s3 = -(1+p_var*(pIndex-1))*(mu)*(s0+dt*s2/2)*(i0+dt*i2/2)/Npop-Vacc1(ii);
    e3 = (1+p_var*(pIndex-1))*(mu)*(s0+(e_var1*v10+e_var2*v20+e_var3*v30)+...
        dt*(s2+(e_var1*v12+e_var2*v22+e_var3*v32))/2)*(i0+dt*i2/2)/Npop - kappa*(e0+dt*e2/2);
    i3 = kappa*(e0+dt*e2/2) - alpha*(i0+dt*i2/2);
    q3 = alpha*(i0+dt*i2/2) - gamma*(q0+dt*q2/2);
    r3 = (1-f)*gamma*(q0+dt*q2/2);
    d3 = f*gamma*(q0+dt*q2/2);
    v13 = Vacc1(ii) - Vacc2(ii) - (1+p_var*(pIndex-1))*(mu)*e_var1*(v10+dt*v12/2)*(i0+dt*i2/2)/Npop;
    v23 = Vacc2(ii) - Vacc3(ii) - (1+p_var*(pIndex-1))*(mu)*e_var2*(v20+dt*v22/2)*(i0+dt*i2/2)/Npop;
    v33 = Vacc3(ii) - (1+p_var*(pIndex-1))*(mu)*e_var3*(v30+dt*v32/2)*(i0+dt*i2/2)/Npop - omega*(v30+dt*v32/2);
    p3 = omega*(v30+dt*v32/2);

    %             if pIndex
    %                 p_var = fun(fit,(x0+(ii-1)/24+dt)/100);
    %             else
    %                 p_var = 0;
    %             end


    s4 = -(1+p_var*(pIndex-1))*(mu)*(s0+dt*s3)*(i0+dt*i3)/Npop-Vacc1(ii);
    e4 = (1+p_var*(pIndex-1))*(mu)*(s0+(e_var1*v10+e_var2*v20+e_var3*v30)+...
        dt*(s3+(e_var1*v13+e_var2*v23+e_var3*v33)))*(i0+dt*i3)/Npop - kappa*(e0+dt*e3);
    i4 = kappa*(e0+dt*e3) - alpha*(i0+dt*i3);
    q4 = alpha*(i0+dt*i3) - gamma*(q0+dt*q3);
    r4 = (1-f)*gamma*(q0+dt*q3);
    d4 = f*gamma*(q0+dt*q3);
    v14 = Vacc1(ii) - Vacc2(ii) - (1+p_var*(pIndex-1))*(mu)*e_var1*(v10+dt*v13)*(i0+dt*i3)/Npop;
    v24 = Vacc2(ii) - Vacc3(ii) - (1+p_var*(pIndex-1))*(mu)*e_var2*(v20+dt*v23)*(i0+dt*i3)/Npop;
    v34 = Vacc3(ii) - (1+p_var*(pIndex-1))*(mu)*e_var3*(v30+dt*v33)*(i0+dt*i3)/Npop - omega*(v30+dt*v33);
    p4 = omega*(v30+dt*v33);

    Y(1,ii+1) = s0+dt*(s1+2*s2+2*s3+s4)/6;
    Y(2,ii+1) = e0+dt*(e1+2*e2+2*e3+e4)/6;
    Y(3,ii+1) = i0+dt*(i1+2*i2+2*i3+i4)/6;
    Y(4,ii+1) = q0+dt*(q1+2*q2+2*q3+q4)/6;
    Y(5,ii+1) = r0+dt*(r1+2*r2+2*r3+r4)/6;
    Y(6,ii+1) = d0+dt*(d1+2*d2+2*d3+d4)/6;
    Y(7,ii+1) = v10+dt*(v11+2*v12+2*v13+v14)/6;
    Y(8,ii+1) = v20+dt*(v21+2*v22+2*v23+v24)/6;
    Y(9,ii+1) = v30+dt*(v31+2*v32+2*v33+v34)/6;
    Y(10,ii+1) = p0+dt*(p1+2*p2+2*p3+p4)/6;
    
end

% Y = round(Y);
%% Write the outputs
S = Y(1,1:N);
E = Y(2,1:N);
I = Y(3,1:N);
Q = Y(4,1:N);
R = Y(5,1:N);
D = Y(6,1:N);
V1 = Y(7,1:N);
V2 = Y(8,1:N);
V3 = Y(9,1:N);
P = Y(10,1:N);


%% Nested functions
   
end


