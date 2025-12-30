%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   EXERCISE 3 - APPLICATIONS OF MATLAB TO THERMAL PROBLEMS
%   Ex 4: Option1 - 4 materials 
%   Abril Bertran
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long
clc;
clear;

%% Numerical Properties

% Mesh
w=1;                        % width of the VC
Lx=1.1;     Ly=0.8;         % Length of each side
Nx=55;      Ny=40;          % Node discretization
dx=Lx/Nx;   dy=Ly/Ny;       % Differentials
ds=dx*w;    dv=dx*dy*w;     

x=zeros(Ny,Nx);
y=zeros(Ny,Nx);

dt=1;                       %time discretization
tend=10000;                  
time=0:dt:tend; time=time.';

Nt=size(time,1);

% Construction of x and y coordinate node matrices
for i=1:Nx
    x(:,i)=0.5*dx+(i-1)*dx;
end

for i=1:Ny
    y(i,:)=0.5*dy+(i-1)*dy;
end
y=flipud(y);

%% Physical Properties

% Material properties
mat1.rho=1500; mat1.cp=750; mat1.lambda=170;
mat2.rho=1600; mat2.cp=770; mat2.lambda=140;
mat3.rho=1900; mat3.cp=810; mat3.lambda=200;
mat4.rho=2500; mat4.cp=930; mat4.lambda=140;

%Construction of material matrices
rho=zeros(Ny,Nx);
cp=zeros(Ny,Nx);
lambda=zeros(Ny,Nx);

for i=1:Nx
    for j=1:Ny
        if x(j,i)<=0.5 && y(j,i)<=0.4
            mat=mat1;
        elseif x(j,i)>0.5 && y(j,i)<=0.7
            mat=mat2;
        elseif x(j,i)<=0.5 && y(j,i)>0.4
            mat=mat3;
        else
            mat=mat4;
        end
        rho(j,i)=mat.rho;
        cp(j,i)=mat.cp;
        lambda(j,i)=mat.lambda;
    end
end


%% Matrices inicialization to save the temperatures
T_all=zeros(Ny,Nx,Nt);              % Matrix temperature to save all data
T_all(:,:,1)=ones(Ny,Nx)*8;         % Temperature at t=0;
T=zeros(Ny,Nx);                     % Matrix temperature

Tp=T_all(:,:,1);                    % Temperature at T^n
T_prev=Tp;                          % Initial estimation for the first iteration


%% Gauss Seidel loop

epsilon=10^-6;

% Time loop to get T^n+1
for t=1:Nt-1
    err=1;
    
    % Iterative loop ensuring convergence
    while err>epsilon
        
        % Calculation to obtain T^n+1    
        for j=Ny:-1:1
            for i=1:Nx
                T(j,i)=getTemperature(i,j,Nx,Ny,lambda,rho,cp,dx,dy,ds,dv,dt,t,T,Tp);
            end
        end
        err=max(max(abs(T-T_prev))); % Error calculation
        T_prev=T;               % Saves the previous results for next iteration
    end
    T_all(:,:,t+1)=T;
    Tp=T;                 % Updates the temperature of T^n for the next timestep iteration

end 


%% Export data

% Obtain the matrix index for the corresponding coordinates
x1=0.65; y1=0.56;
x2=0.74; y2=0.72;

[~, idx_x1] = min(abs(x(1,:) - x1));
[~, idx_y1] = min(abs(y(:,1) - y1));
[~, idx_x2] = min(abs(x(1,:) - x2));
[~, idx_y2] = min(abs(y(:,1) - y2));

% Obtain temperature vectors
T1=T_all(idx_y1,idx_x1,:); T1=T1(:);
T2=T_all(idx_y2,idx_x2,:); T2=T2(:);


% Arrange and export data
data=[time,T1,T2];
dlmwrite('temperaturas.txt', data, 'delimiter', '\t');


%% plot the isotherms for diferent timesteps
t1=2000;
t2= 3000;
t3= 4000;
t4= 5000;

[~, idx_t1] = min(abs(time - t1));
[~, idx_t2] = min(abs(time - t2));
[~, idx_t3] = min(abs(time - t3));
[~, idx_t4] = min(abs(time - t4));

figure
subplot(2, 2, 1);
contourf(x,y,T_all(:,:,idx_t1),'ShowText','on');
set(gca, 'YDir', 'normal');
axis equal;
title('t=2000s');

subplot(2, 2, 2);
contourf(x,y,T_all(:,:,idx_t2),'ShowText','on');
set(gca, 'YDir', 'normal');
axis equal;
title('t=3000s')

subplot(2, 2, 3);
contourf(x,y,T_all(:,:,idx_t3),'ShowText','on');
set(gca, 'YDir', 'normal');
axis equal;
title('t=4000s')

subplot(2, 2, 4);
contourf(x,y,T_all(:,:,idx_t4),'ShowText','on');
set(gca, 'YDir', 'normal');
title('t=5000s')
axis equal;

%% Video Animation

if exist('temperatura_video.mp4','file')
    delete('temperatura_video.mp4');
end

v = VideoWriter('temperatura_video.mp4', 'MPEG-4');
v.FrameRate = 40;  
open(v);

figure('Name','Animation of temprature along time')
ax = axes;
colorbar;
xlabel('x [m]');
ylabel('y [m]');

for t_idx = 2:32:Nt
    contourf(ax, x, y, T_all(:,:,t_idx),'ShowText','on');
    title(['Time: ', num2str(time(t_idx)), ' s']);
    drawnow;
    frame = getframe(gcf); 
    writeVideo(v, frame);
end

close(v);