function [phi_new]=getPhi(rho,gamma,dt,u,v,dx,ds,dv,Nx,Ny,alpha,x,epsilon)
    % Phi inicialization
    
    phi_new= zeros(Ny,Nx);              % Inicialization of the matrix to obtain phi^n+1
    phi_p= zeros (Ny,Nx);                % Initial values of phi at t=0
    
    % Iterative loop ensuring convergencefor t^n
    err=1;
    while err>epsilon
        % Calculation to obtain phi^n+1    
        for j=Ny:-1:1
            for i=1:Nx
                [conv,diff]=getCoeff(i,j,Nx,Ny,gamma,rho,alpha,x,dx,ds,phi_p,u,v);
                phi_new(j,i)=phi_p(j,i)+dt/(rho*dv)*(-conv+diff);
            end
        end
        err=max(max(abs(phi_new-phi_p))); % Error calculation
        phi_p=phi_new;                    % Saves the previous results for next iteration
    end

end