%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   EXERCISE 1 - APPLICATIONS OF MATLAB TO THERMAL PROBLEMS
%   1D transient conduction case - optional part
%   Abril Bertran
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long
clc;
clear;



%% Numerical Temperatures 
vec_dt=linspace(0.001,10,60);     %vector to iterate diferent timesteps
vec_dt1=linspace(0.001,0.4,25);   %vector to iterate diferent timesteps for explicit method
x=0.75;
tend=600;

% Calc of the analytical solution
for i=1:size(vec_dt,2)
    [T_analytical,~]=analytical(vec_dt(i),tend,x);
end

% Explicit beta=0
[Tnum0]=solver(10^-4,0,0.75,tend);

err_dt1=zeros(size(vec_dt1,2),1);
err_dx1=ones(size(vec_dt1,2),1)*abs(Tnum0-T_analytical(end));

for i=1:size(vec_dt1,2)
    [Tnum]=solver(vec_dt1(i),0,x,tend);    
    err_dt1(i)=abs(Tnum0-Tnum);
end


% Implicit beta= 1
[Tnum0]=solver(10^-4,1,0.75,tend);

err_dt2=zeros(size(vec_dt,2),1);
err_dx2=ones(size(vec_dt,2),1)*abs(Tnum0-T_analytical(end));

for i=1:size(vec_dt,2)
    [Tnum]=solver(vec_dt(i),1,x,tend);   
    err_dt2(i)=abs(Tnum0-Tnum);
end

% Crank-Nickolson beta= 0.5
[Tnum0]=solver(10^-4,0.5,0.75,tend);

err_dt3=zeros(size(vec_dt,2),1);
err_dx3=ones(size(vec_dt,2),1)*abs(Tnum0-T_analytical(end));

for i=1:size(vec_dt,2)
    [Tnum]=solver(vec_dt(i),0.5,x,tend); 
    err_dt3(i)=abs(Tnum0-Tnum);
end


%% Plots and Solutions

figure 

loglog(vec_dt1, err_dx1+err_dt1);
hold on
loglog(vec_dt, err_dx2+err_dt2);
loglog(vec_dt, err_dx3+err_dt3);

xlabel('\Delta t');
ylabel('|T-T_{an}|');
title('Discretization Error');
legend('Explicit (\beta=0)','Implicit (\beta=1)','Crank-Nicolson (\beta=0.5)');


hold off
figure

loglog(vec_dt1, err_dt1);
hold on
loglog(vec_dt, err_dt2);
loglog(vec_dt, err_dt3);

xlabel('\Delta t');
ylabel('|T-T_{num}|');
title('Temporal discretization error');
legend('Explicit (\beta=0)','Implicit (\beta=1)','Crank-Nicolson (\beta=0.5)')


[~, idx] = min(abs(vec_dt - 0.1));
disp('error Exp beta=0');
disp(err_dt1(idx));
disp('error Imp beta=1');
disp(err_dt2(idx));
disp('error C-N beta=0.5');
disp(err_dt3(idx));