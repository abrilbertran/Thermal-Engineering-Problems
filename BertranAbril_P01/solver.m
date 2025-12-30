%% FUNCTION GAUSS-SEIDEL
function [Tnum]=solver(dt,beta,x,t)
    L= 1;                     %m
    lambda= 400;              %W/(m*K)
    rho= 8960;                %kg/m^3
    cp= 380;                  %J/(kg*K)
    alpha=lambda/(rho*cp);

    % INPUT NUMERICAL PARAMETERS

    N=100;                                      %num of nodes                            
    dx= L/N;                                    %step x
    vec_x=dx/2:dx:L-dx/2;;                      %discretization of the bar
    vec_t=0:dt:800;                            %vector of time steps
    T= zeros(size(vec_x,2),size(vec_t,2));      %matrix of temperatures
                   
    gamma= alpha*dt/dx^2;

    % BOUNDARY CONDITIONS
    T0=30;        %ºC
    Tl=100;       %ºC
    Tr=20;        %ºC

    % DATA INITIALIZATION
    T(:,1)=T0;      %All the nodes at the initial instant is T0
    T(1,:)=Tl;      %Boundary condition of the first node at all instant =Tl
    T(end,:)=Tr;    %Boundary condition of the last node at all instant =Tr


    % GAUSS - SEIDEL
    tol=1e-13;

    err= zeros(size(vec_t,2),1);

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
        end
    end

    [~, idx1] = min(abs(vec_x - x));
    [~, idx2] = min(abs(vec_t - t));
    Tnum=T(idx1,idx2);

end
