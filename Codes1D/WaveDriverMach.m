% Driver script for solving the 1D advection equations
close all

Globals1D;

% Order of polymomials used for approximation 
N = 6;
k=20;

% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(-pi,pi,k);

% Initialize solver and construct grid and metric
StartUp1D;
% vmapP(1) = vmapM(end);
% vmapP(end) = vmapM(1);

% Set initial conditions
u = cos(x);
p= sin(x);
M=ones(size(x));

FinalTime = 1.5;
[u,p] = Wave1Dmach(u,p,FinalTime,M);
plot(x(:),u(:)); 
hold on
ua= cos(x-FinalTime);
plot(x(:),ua(:))

pa= sin(x-FinalTime);
figure
plot(x(:),p(:)); 
hold on
plot(x(:),pa(:)); 

err = ua - u; % compute point-wise error
M = inv(V*V'); % mass matrix
errL2 = zeros(K,1);
for k = 1 : K
    errL2(k) = err(:,k)'*diag(J(:,k))*M*err(:,k);
end
errL2 = sqrt(sum(errL2)) % Global L^2-norm of error

error=zeros(3,7);
i=1;
for N=[1,2,3,4]
    j=1;
    for kk=2.^(4:10)
        % Generate simple mesh
        [Nv, VX, K, EToV] = MeshGen1D(-pi,pi,kk);
        
        % Initialize solver and construct grid and metric
        StartUp1D;
        
        % Set initial conditions
        u = cos(x);
        p= sin(x);
        M=ones(size(x));
        
        FinalTime = 1.5;
        [u,p] = Wave1Dmach(u,p,FinalTime,M);
        
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
    i=i+1
end

%%
N=[1,2,3,4];

figure
for i=1:4
    loglog(2.^(4:10),error(i,:),'-*','LineWidth',2)
    hold on
end

for i=1:4
    loglog(2.^(4:10),1*(2.^(4:10)).^(-(N(i)+1)),'--','LineWidth',2)
end


xlabel('$K$', 'Interpreter', 'Latex', 'FontSize', 15)
ylabel('$\| \varepsilon_h(T)\|$', 'Interpreter', 'Latex','FontSize', 15)
legend('$N=1$','$N=2$','$N=3$','$N=4$','$\mathcal{O}(h^{N+1})$','Interpreter', 'Latex','FontSize', 15,'Location','southwest')

