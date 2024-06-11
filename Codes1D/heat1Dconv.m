clear all
close all
clc

error=zeros(4,4);

i=1;
for N=[1,2,3,4]
    j=1;
    for kk=2.^(4:8)
        % Generate simple mesh
        [Nv, VX, K, EToV] = MeshGen1D(-1,1,kk);
        
        % Initialize solver and construct grid and metric
        StartUp1D;
        
        % Set initial conditions
        u = sin(x*pi);

        % Solve Problem
        FinalTime = 0.2;
        [u,time] = Heat1D(u,FinalTime);

        ua=exp(-pi^2*FinalTime)*sin(pi*x);
        
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
for i=1:4
    loglog(2.^(4:8),error(i,:),'-*','linewidth',2)
    hold on
end
loglog(2.^(4:8),0.1*(2.^(4:8)).^(-1),'--','linewidth',2)
loglog(2.^(4:8),0.1*(2.^(4:8)).^(-3),'--','linewidth',2)
loglog(2.^(4:8),0.1*(2.^(4:8)).^(-5),'--','linewidth',2)
legend('$N=1$','$N=2$','$N=3$','$N=4$','$\mathcal{O}(h^{1})$','$\mathcal{O}(h^3)$','$\mathcal{O}(h^{5})$','location','southwest','Interpreter', 'Latex', 'FontSize', 15)
xlabel('$K$', 'Interpreter', 'Latex', 'FontSize', 15)
ylabel('$\| \varepsilon_h(T)\|$', 'Interpreter', 'Latex','FontSize', 15)
