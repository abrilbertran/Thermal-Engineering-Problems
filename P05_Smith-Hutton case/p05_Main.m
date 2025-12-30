%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   EXERCISE 5 - APPLICATIONS OF MATLAB TO THERMAL PROBLEMS
%   Smith Hutton Case
%   Abril Bertran
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long
clc;
clear;

%% PHYSICAL PROPERTIES
L=2;
H=1;
gamma=1;
rho=10;
alpha=10;

%% NUMERICAL PROPERTIES

%Time discretization
dt=0.001;           

w=1;                      % width of the VC
Nx=100;    Ny=50;         % Node discretization
dx=L/Nx;   dy=H/Ny;       % Differentials
ds=dx*w;   dv=dx*dy*w;

% Construction of x and y coordinate node matrices
x=zeros(Ny,Nx);
y=zeros(Ny,Nx);

for i=1:Nx
    x(:,i)=-1+0.5*dx+(i-1)*dx;
end

for i=1:Ny
    y(i,:)=0.5*dy+(i-1)*dy;
end
y=flipud(y);


%% Velocity field
u=2*y.*(1-x.^2);
v=-2*x.*(1-y.^2);

% Plot of the velocity
%quiver(x,y,u,v,2);
%axis equal;


%% Phi inicialization

phi_new= zeros(Ny,Nx);              % Inicialization of the matrix to obtain phi^n+1
phi_p= zeros (Ny,Nx);               % Initial values of phi at t=0

%% Phi calculation

epsilon=10^-6;
err=1;
iter=0;

[phi_1]=getPhi(10,gamma,dt,u,v,dx,ds,dv,Nx,Ny,alpha,x,epsilon);
[phi_2]=getPhi(1000,gamma,dt,u,v,dx,ds,dv,Nx,Ny,alpha,x,epsilon);
[phi_3]=getPhi(10^5,gamma,10^-5,u,v,dx,ds,dv,Nx,Ny,alpha,x,epsilon);



%% Data visualization
T = table(x(end,:).', phi_1(end,:).',phi_2(end,:).',phi_3(end,:).', 'VariableNames', {'x', 'rho10','rho1000','rho10e6'});
T=T(T.x>0,:);

disp(T);

%% plot the color map

figure
surf(x,y,phi_1);
xlim([-1 1]);
ylim([0 2]);
axis equal
title('\rho / \Gamma=10');
shading interp;       % suaviza el gráfico
view(2);              
colorbar;

figure
surf(x,y,phi_2);
xlim([-1 1]);
ylim([0 2]);
axis equal 
title('\rho / \Gamma=1000');
shading interp;       
view(2);              
colorbar;

figure
surf(x,y,phi_3);
xlim([-1 1]);
ylim([0 2]);
axis equal
title('\rho / \Gamma=1000000');
shading interp;       
view(2);              
colorbar;

%% Plot the outlet profile
figure
subplot(1,3,1);
plot(x(end,:), phi_1(end,:),'LineWidth', 1.5);
xlabel('x'); 
ylabel('\phi');
ylim ([0 2.1]);
title('\rho / \Gamma=10');

subplot(1,3,2);
plot(x(end,:), phi_2(end,:),'LineWidth', 1.5);
xlabel('x'); 
ylabel('\phi');
ylim ([0 2.1]);
title('\rho / \Gamma=1000');

subplot(1,3,3);
plot(x(end,:), phi_3(end,:),'LineWidth', 1.5);
xlabel('x'); 
ylabel('\phi');
ylim ([0 2.1]);
title('\rho / \Gamma=1000000');




