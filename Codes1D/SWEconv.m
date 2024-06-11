clear all
close all
clc



error=zeros(3,3);

i=1;
for N=[1,2,3]
    j=1;
    for kk=2.^(4:8)
        Globals1D;
        A=1;
        l=2*pi;
        T=2*pi;
        omega=2*pi/T;
        k=2*pi/l;
        g=9.81;
        h0=omega^2/(k^2*g);

        % Generate simple mesh
       [Nv, VX, K, EToV] = MeshGen1D(0,2*pi,kk);
        
        % Initialize solver and construct grid and metric
        StartUp1D;
        vmapP(1) = vmapM(end);
        vmapP(end) = vmapM(1);
        
        eta = A*cos(omega*0-k*x);
        u= omega/(k*h0)*eta;

        FinalTime = 5;
        [eta,u] = SWE1D(eta,u,FinalTime);

        etaa= A*cos(omega*FinalTime-k*x);
        ua=omega/(k*h0)*etaa;
        
        err = ua - u; % compute point-wise error
        M = inv(V*V'); % mass matrix
        errL2 = zeros(K,1);
        for k = 1 : K
            errL2(k) = err(:,k)'*diag(J(:,k))*M*err(:,k);
        end
        error(i,j) = sqrt(sum(errL2)); % Global L^2-norm of error
        j=j+1;
        
    end
    i=i+1
end

%%
figure
for i=1:3
    loglog(2.^(4:8),error(i,:),'-*','linewidth',2)
    hold on
end
loglog(2.^(4:8),5*(2.^(4:8)).^(-1),'--','linewidth',2)
loglog(2.^(4:8),10*(2.^(4:8)).^(-3),'--','linewidth',2)
% loglog(2.^(4:8),15*(2.^(4:8)).^(-3),'--','linewidth',2)
legend('$N=1$','$N=2$','$N=3$','$\mathcal{O}(h^{2})$','$\mathcal{O}(h^3)$','$\mathcal{O}(h^{4})$','location','southwest','Interpreter', 'Latex', 'FontSize', 15)
xlabel('$K$', 'Interpreter', 'Latex', 'FontSize', 15)
ylabel('$\| \varepsilon_h(T)\|$', 'Interpreter', 'Latex','FontSize', 15)
