function [T_analytical,vec_t]=analytical(dt,tend,x)

    % INPUT PHYSICAL PARAMETERS
    L= 1;                     %m
    lambda= 400;              %W/(m*K)
    rho= 8960;                %kg/m^3
    cp= 380;                  %J/(kg*K)
    alpha=lambda/(rho*cp);
    
    % BOUNDARY CONDITIONS
    T0=30;        %ºC
    Tl=100;       %ºC
    Tr=20;        %ºC
    
    
    % ANALYTICAL SOLUTION
    
    % Calculation of the temperature at given x along the time
    vec_t=0:dt:tend;
    x_adim=x/L;
    t_adim= vec_t.'*(alpha/L^2);
    u0=(T0-Tl)/(Tr-Tl);
    u1=1;
    
    u=zeros(size(vec_t,2),1);
    u=1*x_adim+u;               %Calculation of the initial term
    
    
    % Sum of impair index
    sum1=zeros(size(vec_t,2),1);
    for i=1:2:100
        sum1=sum1+sin(i*pi*x_adim)/i*exp(-i^2*pi^2.*t_adim);
    end
    
    %Sum of all index
    sum2=zeros(size(vec_t,2),1);
    for i=1:100
        sum2=sum2+(-1)^i*sin(i*pi*x_adim)/i*exp(-i^2*pi^2.*t_adim);
    end
    
    u=u+4*u0/pi*sum1+2*u1/pi*sum2;
    T_analytical=Tl+u*(Tr-Tl);

end