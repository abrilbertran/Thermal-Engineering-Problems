%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   EXERCISE 1 - APPLICATIONS OF MATLAB TO THERMAL PROBLEMS
%   1D transient conduction case - mandatory part
%   Abril Bertran
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long
clc;
clear;

%% INPUT PHYSICAL PARAMETERS

L= 1;                     %m
lambda= 400;              %W/(m*K)
rho= 8960;                %kg/m^3
cp= 380;                  %J/(kg*K)
alpha=lambda/(rho*cp);

%% INPUT NUMERICAL PARAMETERS

N=100;                                      %num of nodes                            
dx= L/N;                                    %step x
vec_x=dx/2:dx:L-dx/2;                      %discretization of the bar

dt= 0.25;                                   %time step
tend=2000;                                  %final time
vec_t=0:dt:tend;                            %vector of time steps
T= zeros(size(vec_x,2),size(vec_t,2));      %matrix of temperatures

beta=0.5;                   
gamma= alpha*dt/dx^2;

%% BOUNDARY CONDITIONS

T0=30;        %ºC
Tl=100;       %ºC
Tr=20;        %ºC

%% DATA INITIALIZATION
T(:,1)=T0;      %All the nodes at the initial instant is T0
T(1,:)=Tl;      %Boundary condition of the first node at all instant =Tl
T(end,:)=Tr;    %Boundary condition of the last node at all instant =Tr


%% GAUSS - SEIDEL
tol=1e-6;
contador= zeros(size(vec_t,2),1);
err= zeros(size(vec_t,2),1);

tic
for n=1:size(T,2)-1
    err(n)=10;

    %Initial estimation
    T(:,n+1)=T(:,n);

    while err(n)>tol
        Told=T(:,n+1);

        for i=2:size(T,1)-1
            Te=T(i+1,n+1);
            Tw=T(i-1,n+1);
    
            %Calcul coeficients
            ap=1+2*gamma*beta;
            ae=gamma*beta;
            aw=gamma*beta;
            bp=T(i,n)+gamma*(1-beta)*(T(i+1,n)-2*T(i,n)+T(i-1,n));

            %Calcul temperatura
            T(i,n+1)=1/ap*(ae*Te+aw*Tw+bp);
        end 

        err(n)= max(abs(Told-T(:,n+1)));
        contador(n)=contador(n)+1;
    end

end
temps=toc;
disp('Temps:');
disp(temps);

%% Plots and Solutions

% 1. Plot for x=0.5
[~, x] = min(abs(vec_x - 0.5)); % Gets the cell index for x=0.5
[T_analytical,~]=analytical(dt,tend,0.5);
figure

plot(vec_t,T(x,:));
hold on
plot(vec_t,T_analytical);
title('Temperatures at x=0.5 along the time')
xlabel('Time [s]');
ylabel('Temperatures [ºC]');
ylim([25,60]);
legend('Numerical Solution','Analytical Solution')

% 2. Temperature profile at t=800s of along x
[~, t] = min(abs(vec_t - 800));   % Gets the cell index for t=800
figure
plot(vec_x,T(:,t));
title('Temperatures at t=800s along the bar')
xlabel('X [m]');
ylabel('Temperatures [ºC]');

% 3. Temperature at x=0.75m and t=600s
[~, x] = min(abs(vec_x - 0.75));
[~, t] = min(abs(vec_t - 600)); 
[T_analytical,~]=analytical(dt,600,0.75);

disp('Numerical Temperature at x=0.5m and t=600s:');
disp(T(x,t));
disp('Analytical Temperature at x=0.5m and t=600s:');
disp(T_analytical(end));