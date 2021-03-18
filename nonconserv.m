function [rho,v,T,k,netmf,err_v,err_rho,err_T] = nonconserv(n,x,dx,gamma,nt,c,tol, mass_tol)

%Intial Profiles
rho = (1 - 0.3146*x); %Density
T = (1 - 0.2314*x); %Temperature
v = (0.1 + (1.09*x)).*(T.^0.5); %Velocity

a = (1 + 2.2*(x-1.5).^2); %Area
% Outer Time Loop
for k = 1:nt
    
    dt = min(abs(c*dx./(realsqrt(T)+v)));
    
    %Copying old variables
    rho_old = rho;
    v_old = v;
    T_old = T;
    
    
    % Predictor Method
    for j = 2:n-1
        dvdx = (v(j+1) - v(j))/dx;
        drhodx = (rho(j+1) - rho(j))/dx;
        dlogadx = (log(a(j+1)) - log(a(j)))/dx;
        dTdx = (T(j+1) - T(j))/dx;
        
        %Continuity Equation
        drhodt_p(j) = -rho(j)*dvdx - rho(j)*v(j)*dlogadx - v(j)*drhodx;
        
        %Momentum Equation
        dvdt_p(j) = -v(j)*dvdx - (1/gamma)*(dTdx + (T(j)/rho(j))*drhodx);
        
        %Energy Equation
        dTdt_p(j) = -v(j)*dTdx - (gamma-1)*T(j)*(dvdx + v(j)*dlogadx);
        
        %Solution Update
        v(j) = v(j) + dvdt_p(j)*dt;
        rho(j) = rho(j) + drhodt_p(j)*dt;
        T(j) = T(j) + dTdt_p(j)*dt;
    end
    
    % Corrector Method
    for j = 2:n-1
        dvdx = (v(j) - v(j-1))/dx;
        drhodx = (rho(j) - rho(j-1))/dx;
        dlogadx = (log(a(j)) - log(a(j-1)))/dx;
        dTdx = (T(j) - T(j-1))/dx;
        
        %Continuity Equation
        drhodt_c(j) = -rho(j)*dvdx - rho(j)*v(j)*dlogadx - v(j)*drhodx;
        
        %Momentum Equation
        dvdt_c(j) = -v(j)*dvdx - (1/gamma)*(dTdx + (T(j)/rho(j))*drhodx);
        
        %Energy Equation
        dTdt_c(j) = -v(j)*dTdx - (gamma-1)*T(j)*(dvdx + v(j)*dlogadx);
    end
    
    %Computing average time derivative
    drhodt = 0.5*(drhodt_p + drhodt_c);
    dvdt = 0.5*(dvdt_p + dvdt_c);
    dTdt = 0.5*(dTdt_p + dTdt_c);
    
    %Final Update
    for j = 2:n-1
        v(j) = v_old(j) + dvdt(j)*dt;
        rho(j) = rho_old(j) + drhodt(j)*dt;
        T(j) = T_old(j) + dTdt(j)*dt;
    end
    
    %Applying Boundary Conditions
    %inlet
    v(1) = 2*v(2) - v(3);
    
    %outlet
    v(n) = 2*v(n-1) - v(n-2);
    rho(n) = 2*rho(n-1) - rho(n-2);
    T(n) = 2*T(n-1) - T(n-2);
    
    %Mass Flow Rate
    inflow = rho(1)*v(1)*a(1);
    outflow = rho(n)*v(n)*a(n);
    netmf(k) = (outflow-inflow);
    
    %Error Tolerance
    err_v(k) = max(v - v_old);
    err_rho(k) = max(rho - rho_old);
    err_T(k) = max(T - T_old);
    if(err_v(k) < tol && err_rho(k) < tol && err_T(k) < tol && netmf(k)<mass_tol)
        break;
    end
end
end