% Driver script for solving the 1D advection equations
close all

Globals1D;

% Order of polymomials used for approximation 

error=zeros(2,7);

i=1;
for N=[1,2]
    j=1;
    for kk=2.^(4:11)
        % Generate simple mesh
        [Nv, VX, K, EToV] = MeshGen1D(-2,2,kk);
        
        % Initialize solver and construct grid and metric
        StartUp1D;
        
        % Set initial conditions
        u= cos(x);
        p= sin(x);
        M=ones(size(x));

        FinalTime = 1;
        [u,p] = Wave1Dmach(u,p,FinalTime,M);
        plot(x(:),u(:)); 
        hold on
        ua= cos(x-FinalTime);
        
        err = ua - u; % compute point-wise error
        M = inv(V*V'); % mass matrix
        errL2 = zeros(K,1);
        for k = 1 : K
            errL2(k) = err(:,k)'*diag(J(:,k))*M*err(:,k);
        end
        error(i,j) = sqrt(sum(errL2)); % Global L^2-norm of error
        j=j+1;
       
    end
    i=i+1;
end

%%
N=[1,2]
figure
for i=1:2
    loglog(2.^(4:11),error(i,:),'-*','LineWidth',2)
    hold on
end

for i=1:3
    loglog(2.^(4:11),1*(2.^(4:11)).^(-(N(i)+1)),'--','LineWidth',2)
end


xlabel('$K$', 'Interpreter', 'Latex', 'FontSize', 15)
ylabel('$\| \varepsilon_h(T)\|$', 'Interpreter', 'Latex','FontSize', 15)
legend('$N=1$','$N=2$','$N=4$','$N=8$','$\mathcal{O}(h^{N+1})$','Interpreter', 'Latex','FontSize', 15,'Location','southwest')

