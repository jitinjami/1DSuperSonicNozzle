clear all
close all
clc

%Inputs
n = 31; %Number of nodes
x = linspace(0,3,n); %Mesh
dx = x(2) - x(1);
gamma = 1.4;

%Time steps
nt = 5000;
c1 = 0.5;
c2 = 0.5;
tol = 1e-6;
mass_tol1 = 1e-3;
mass_tol2 = 1e-2;

[rho1,v1,T1,total_time1,netmf1,err_v1,err_rho1,err_T1] = nonconserv(n,x,dx,gamma,nt,c1,tol,mass_tol1);

[rho2,v2,T2,total_time2,netmf2,err_v2,err_rho2,err_T2] = conserv(n,x,dx,gamma,nt,c2,tol,mass_tol2);

%Velocity plots
figure(1)
plot(x,v1,'b',x,v2,'r')
xlabel('Non-dimentional Distance')
ylabel('Non-dimentional Velocity')
legend('Non Conservative','Conservative')
grid on


%Density plots
figure(2)
plot(x,rho1,'b',x,rho2,'r')
xlabel('Non-dimentional Distance')
ylabel('Non-dimentional Density')
legend('Non Conservative','Conservative')
grid on

%Temperature plots
figure(3)
plot(x,T1,'b',x,T2,'r')
xlabel('Non-dimentional Distance')
ylabel('Non-dimentional Temperature')
legend('Non Conservative','Conservative')
grid on

%Net Mass Flow Rate Plot
figure(4)
subplot(2,1,1)
plot(netmf1)
xlabel('Time Step')
ylabel('Net Mass Flow Rate')
title('Non Conservative Form')
grid on
subplot(2,1,2)
plot(netmf2)
xlabel('Time Step')
ylabel('Net Mass Flow Rate')
title('Conservative Form')
grid on

%Error in Velocity with Time
figure(5)
subplot(2,1,1)
plot(err_v1)
xlabel('Time Step')
ylabel('Error in Velocity')
title('Non Conservative Form')
grid on
subplot(2,1,2)
plot(err_v2)
xlabel('Time Step')
ylabel('Error in Velocity')
title('Conservative Form')
grid on

%Error in density with Time
figure(6)
subplot(2,1,1)
plot(err_rho1)
axis([0 inf -0.01 0.02])
xlabel('Time Step')
ylabel('Error in Density')
title('Non Conservative Form')
grid on
subplot(2,1,2)
plot(err_rho2)
axis([0 inf -0.005 0.015])
xlabel('Time Step')
ylabel('Error in Density')
title('Conservative Form')
grid on

%Error in Temperature with Time
figure(7)
subplot(2,1,1)
plot(err_T1)
axis([0 inf -0.005 0.01])
xlabel('Time Step')
ylabel('Error in Temperature')
title('Non Conservative Form')
grid on
subplot(2,1,2)
plot(err_T2)
axis([0 inf -0.01 0.02])
xlabel('Time Step')
ylabel('Error in Temperature')
title('Conservative Form')
grid on