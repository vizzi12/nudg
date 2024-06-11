clear all
close all
clc

error=zeros(3,3);

i=1;
for N=[1,2,3]
    j=1;
    for kk=2.^(4:7)
        % Generate simple mesh
        [Nv, VX, K, EToV] = MeshGen1D(-1,1,kk);
        
        % Initialize solver and construct grid and metric
        StartUp1D;
       

        % Set initial conditions
        u = cos(x*pi);
        FinalTime = 0.0025;
        [u,time] = Dispersive1D(u,FinalTime);

        ua=cos(pi^3*FinalTime+pi*x);
        
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
    loglog(2.^(4:7),error(i,:),'-*','linewidth',2)
    hold on
end
loglog(2.^(4:7),0.2*(2.^(4:7)).^(-1),'--','linewidth',2)
loglog(2.^(4:7),0.2*(2.^(4:7)).^(-3),'--','linewidth',2)
legend('$N=1$','$N=2$','$N=3$','$\mathcal{O}(h^1)$','$\mathcal{O}(h^3)$','location','southwest','Interpreter', 'Latex', 'FontSize', 15)
xlabel('$K$', 'Interpreter', 'Latex', 'FontSize', 15)
ylabel('$\| \varepsilon_h(T)\|$', 'Interpreter', 'Latex','FontSize', 15)
