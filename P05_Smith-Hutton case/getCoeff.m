function [conv,diff]=getCoeff(i,j,Nx,Ny,gamma,rho,alpha,x,dx,ds,phi_p,u,v)
    
    % Velocities
    if i>1,  uw=0.5*(u(j,i)+u(j,i-1)); else, uw=u(j,i); end
    if i<Nx, ue=0.5*(u(j,i)+u(j,i+1)); else, ue=u(j,i); end
    if j>1,  vn=0.5*(v(j,i)+v(j-1,i)); else, vn=v(j,i); end
    if j<Ny, vs=0.5*(v(j,i)+v(j+1,i)); else, vs=v(j,i); end

    
    % Calculation of phis by CDS
    if i>1,  phiW=0.5*(phi_p(j,i)+phi_p(j,i-1)); else, phiW=1-tanh(alpha); end
    if i<Nx, phiE=0.5*(phi_p(j,i)+phi_p(j,i+1)); else, phiE=1-tanh(alpha); end
    if j>1,  phiN=0.5*(phi_p(j,i)+phi_p(j-1,i)); else, phiN=1-tanh(alpha); end
    if j<Ny, phiS=0.5*(phi_p(j,i)+phi_p(j+1,i)); end


    % Bottom wall
    if j==Ny
        if x(j,i) <=0  %inlet
            phiS= 1+tanh((2*x(j,i)+1)*alpha);
            conv=rho*ds*(ue*phiE-uw*phiW+vn*phiN-vs*phiS);
            diff=gamma*ds/dx*(phiE+phiW+phiN+phiS-4*phi_p(j,i));

        else %outlet
            phiS=phi_p(j,i);
            conv=rho*ds*(ue*phiE-uw*phiW+vn*phiN-vs*phiS);
            diff=gamma*ds/dx*(phiE+phiW+phiN-3*phi_p(j,i));
        end
        return  %finalizes the execution of the function, useful for this specific case
    end


    % General case
    conv=rho*ds*(ue*phiE-uw*phiW+vn*phiN-vs*phiS);
    diff=gamma*ds/dx*(phiE+phiW+phiN+phiS-4*phi_p(j,i));


end