function [rho,v,T,k,netmf,err_v,err_rho,err_T] = conserv(n,x,dx,gamma,nt,c,tol,mass_tol)

%Intial Profiles
rho = 1 - 0.3146*x; %Density
T = 1 - 0.2312*x; %Temperature
a = 1 + 2.2*(x-1.5).^2; %Area
v = 0.59./(rho.*a); %Velocity

%Derivatives
U1 = rho.*a;
U2 = rho.*a.*v;
U3 = rho.*(T/(gamma-1) + (gamma/2)*v.^2).*a;


% Outer Time Loop
for k = 1:nt
    
    % storing old values
    U1_old = U1;
    U2_old = U2;
    U3_old = U3;
    
    rho_old = U1_old./a;
    v_old = U2_old./(rho.*a);
    T_old = (gamma-1)*((U3_old./U1_old) - (gamma/2)*v.^2);
    
    F1 = U2;
    F2 = (U2.^2)./U1 + ((gamma-1)/gamma)*(U3 - (gamma/2)*(U2.^2)./U1);
    F3 = (gamma*U2.*U3)./U1 - ((gamma*(gamma-1)/2)*(U2.^3)./U1.^2);
    
    for j = 2:n-1
        J2(j) = (1/gamma)*rho(j)*T(j)*(a(j+1) - a(j))/dx;
    end
    
    for i = 1:n
        deltat(i)=(c*(dx/((T(i)^0.5+v(i)))));
    end
    dt = min(deltat);
    
    % Predictor Method
    for j = 2:n-1
        %Continuity Equation
        dU1dt_p(j) = -(F1(j+1)-F1(j))/dx;
        
        %Momentum Equation
        dU2dt_p(j) = -(F2(j+1)-F2(j))/dx + J2(j);
        
        %Energy Equation
        dU3dt_p(j) = -(F3(j+1)-F3(j))/dx;
        
        %Solution Update
        U1(j) = U1(j) + dU1dt_p(j)*dt;
        U2(j) = U2(j) + dU2dt_p(j)*dt;
        U3(j) = U3(j) + dU3dt_p(j)*dt;
    end
    
    %Calculating Primitive variables
    rho = U1./a;
    v = U2./(rho.*a);
    T = (gamma-1)*((U3./U1) - (gamma/2)*v.^2);
    
    F1 = U2;
    F2 = (U2.^2)./U1 + ((gamma-1)/gamma)*(U3 - (gamma/2)*(U2.^2)./U1);
    F3 = (gamma*U2.*U3)./U1 - ((gamma*(gamma-1)/2)*(U2.^3)./U1.^2);
    
    for j = 2:n-1
        J2(j) = (1/gamma)*rho(j)*T(j)*(a(j) - a(j-1))/dx;
    end
    
    % Corrector Method
    for j = 2:n-1
        %Continuity Equation
        dU1dt_c(j) = -(F1(j)-F1(j-1))/dx;
        
        %Momentum Equation
        dU2dt_c(j) = -(F2(j)-F2(j-1))/dx + J2(j);
        
        %Energy Equation
        dU3dt_c(j) = -(F3(j)-F3(j-1))/dx;
    end
    
    %Computing average time derivative
    dU1dt = 0.5*(dU1dt_p + dU1dt_c);
    dU2dt = 0.5*(dU2dt_p + dU2dt_c);
    dU3dt = 0.5*(dU3dt_p + dU3dt_c);
    
    %Final Update
    for j = 2:n-1
        U1(j) = U1_old(j) + dU1dt(j)*dt;
        U2(j) = U2_old(j) + dU2dt(j)*dt;
        U3(j) = U3_old(j) + dU3dt(j)*dt;
    end
    
    %Applying Boundary Conditions
    %inlet
    U2(1) = 2*U2(2) - U2(3);
    
    %outlet
    U1(n) = 2*U1(n-1) - U1(n-2);
    U2(n) = 2*U2(n-1) - U2(n-2);
    U3(n) = 2*U3(n-1) - U3(n-2);
    
    %Recalculating Primitive variables
    rho = U1./a;
    v = U2./(rho.*a);
    T = (gamma-1)*((U3./U1) - (gamma/2)*v.^2);
    
    %Mass Flow Rate
    inflow = rho(1)*v(1)*a(1);
    outflow = rho(n)*v(n)*a(n);
    netmf(k) = outflow-inflow;
    
    %Error Tolerance
    err_v(k) = max(v - v_old);
    err_rho(k) = max(rho - rho_old);
    err_T(k) = max(T - T_old);
    if(err_v(k) < tol && err_rho(k) < tol && err_T(k) < tol && netmf(k)<mass_tol)
        break;
    end
end
end