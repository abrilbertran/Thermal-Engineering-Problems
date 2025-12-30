function [Tnew]=getTemperature(i,j,Nx,Ny,lambda,rho,cp,dx,dy,ds,dv,dt,t,T,Tp)

    %% PARAMETERS DEFINITION
    % Boundary conditions
    alpha=9;
    Tbottom=23;
    Tright=8+0.005*t*dt;
    Tleft=33;
    Qtop=60;

    % Parameters for actual node P
    rhoP=rho(j,i);
    cpP=cp(j,i);

    % Calculation of lambdas and harmonic mean
    if i>1, lambW=1/(1/lambda(j,i)+1/lambda(j,i-1)); else, lambW=lambda(j,i); end
    if i<Nx, lambE=1/(1/lambda(j,i)+1/lambda(j,i+1)); else, lambE=lambda(j,i); end
    if j>1, lambN=1/(1/lambda(j,i)+1/lambda(j-1,i)); else, lambN=lambda(j,i); end
    if j<Ny, lambS=1/(1/lambda(j,i)+1/lambda(j+1,i)); else, lambS=lambda(j,i); end
    


    %% COEFFICIENTS GENERAL CASE ( valid for all internal nodes)
    ae=lambE*ds/dx;
    aw=lambW*ds/dx;
    an=lambN*ds/dy;
    as=lambS*ds/dy;
    ap=rhoP*cpP*dv/dt+(ae+aw+an+as);
    b=(rhoP*cpP*dv/dt)*Tp(j,i);


    %% CASE STUDIES

    % Left bottom corner  --> aw=as=0
    if i==1 && j==Ny
        ap=ap-aw+alpha*ds;
        b=b+alpha*Tleft*ds+lambS*ds/dy*Tbottom;
        Tnew=1/ap*(ae*T(j,i+1)+an*T(j-1,i)+b);

    % Right bottom corner --> ae=as=0
    elseif i==Nx && j==Ny
        b=b+lambE*ds/dx*Tright+lambS*ds/dy*Tbottom;
        Tnew=1/ap*(aw*T(j,i-1)+an*T(j-1,i)+b);

    % Left upper corner --> aw=an=0
    elseif i==1 && j==1
        ap=ap-aw-an+alpha*ds;
        b=b+Qtop*dx+alpha*Tleft*ds;
        Tnew=1/ap*(ae*T(j,i+1)+as*T(j+1,i)+b);

    % Right upper corner --> ae=an=0
    elseif i==Nx && j==1
        ap=ap-an;
        b=b+Qtop*dx+lambE*ds/dx*Tright;
        Tnew=1/ap*(aw*T(j,i-1)+as*T(j+1,i)+b);

    % Left side nodes --> aw=0
    elseif i==1
        ap=ap-aw+alpha*ds;
        b=b+alpha*Tleft*ds;
        Tnew=1/ap*(ae*T(j,i+1)+an*T(j-1,i)+as*T(j+1,i)+b);

    % Right side nodes --> ae=0
    elseif i==Nx
        b=b+lambE*ds/dx*Tright;
        Tnew=1/ap*(aw*T(j,i-1)+an*T(j-1,i)+as*T(j+1,i)+b);

    % Bottom side nodes --> as=0
    elseif j==Ny
        b=b+lambS*ds/dy*Tbottom;
        Tnew=1/ap*(ae*T(j,i+1)+aw*T(j,i-1)+an*T(j-1,i)+b);

    % Upper side nodes --> an=0
    elseif j==1
        ap=ap-an;
        b=b+Qtop*dx;
        Tnew=1/ap*(ae*T(j,i+1)+aw*T(j,i-1)+as*T(j+1,i)+b);

    % Internal nodes
    else  
        Tnew=1/ap*(ae*T(j,i+1)+aw*T(j,i-1)+an*T(j-1,i)+as*T(j+1,i)+b);
    end

end