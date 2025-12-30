%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   EXERCISE 2 - APPLICATIONS OF MATLAB TO THERMAL PROBLEMS
%   MATLAB Linear Solvers
%   Abril Bertran
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long
clc;
clear;

%% Physical propperties
T0 = 30;    % Initial temperature [ºC]
Tl = 100;   % Left temperature [ºC]
Tr = 20;    % Right temperature [ºC]
L = 1;      % Bar length [m]
tfin = 600; % End time [s]
alpha = 400/(8960*380); % thermal diffusivity [m2/s]

%% Numerical propperties
N = 500;    % Number of nodes
dt = 0.1;   % Timestep
dx = L/N;



%% Dimension set of matrix A ,vector b and vector T (started at T0)
A = zeros(N);
b = zeros(N,1);
T = ones(N,1)*T0;

%% Dimensio set for vectors P and R
P=zeros(N,1);
R=zeros(N,1);

%% Other calculations
x=dx/2:dx:L-dx/2;
gamma= alpha*dt/dx^2;

%% Construction of A matrix
for i=2:(N)
    for j=2:(N)
        if i==j
            A(i,j)=1+2*gamma;
            A(i-1,j)=-gamma;
            A(i,j-1)=-gamma;
        end
    end
end

A(1,1)=1+3*gamma;
A(end,end)=1+3*gamma;

%% Construction of S matrix
S=sparse(A);

%% Inverted A matrix

invA=inv(A);


%% Temporal Iteration
t = 0;

tic
while t <= tfin
    
    for i=1:length(b)
        if i==1
            b(i)=T(i)+2*gamma*Tl;
        elseif i==length(b)
            b(i)=T(i)+2*gamma*Tr;
        else
            b(i)=T(i);
        end
    end

    % A) MLDIVIDE
    %T_result=mldivide(A,b);

    % B) LINSOLVE
    %opts.POSDEF=true;
    %opts.SYM=true;

    %T_result= linsolve(A,b,opts);

    % C)INVERT MATRIX
    %T_result=inv(A)*b;

    % D) INVERT MATRIX BEFORE TIMELOOP
    %T_result=invA*b;

    % E) SPARSE MATRIX
    %T_result=mldivide(S,b);

    %F) TDMA
    P(1)=A(1,2)/A(1,1);
    R(1)=b(1)/A(1,1);
    for i=2:(N-1)
        P(i)=A(i,i+1)/(A(i,i)-A(i,i-1)*P(i-1));
        R(i)=(b(i)-A(i,i-1)*R(i-1))/(A(i,i)-A(i,i-1)*P(i-1));
    end

    R(N)=(b(N)-A(N,N-1)*R(N-1))/(A(N,N)-A(N,N-1)*P(N-1));
    for i=N:-1:1
        if i==N
            T_result(i)=R(i);
        else
            T_result(i)=R(i)-P(i)*T_result(i+1);
        end
    end

    T=T_result;
    t = t + dt;
end


temps=toc;









