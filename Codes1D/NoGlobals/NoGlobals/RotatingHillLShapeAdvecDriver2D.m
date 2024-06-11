% % Driver script for solving the 1D advection equations
% clear all
% clc
% close all
% 
% % Order of polymomials used for approximation 
% G.N = 5;
% h = 0.2;
% 
% figure1 = figure;
% fd=@(p) ddiff(drectangle(p,-1,1,-1,1),drectangle(p,0,1,0,1));
% [p,G.EToV]=distmesh2d(fd,@huniform,h,[-1,-1;1,1],[-1,-1;-1,1;1,-1;0,1;1,0;0,0]);
% close(figure1)
% 
% G.VX = p(:,1)';
% G.VY = p(:,2)';
% 
% G.EToV = Reorder(G.EToV, G.VX, G.VY);
% 
% G.K = height(G.EToV);
% G.Nv = length(G.VX);
% 
% % Initialize solver and construct grid and metric
% G = StartUp2D(G);
% 
% % Set initial conditions
% cx = 2*pi*G.y;
% cy = -2*pi*G.x;
% lam = 1/8;
% u = Advec2Dsol(G.x, G.y, 0, lam);
% 
% % Solve Problem
% FinalTime = 1;
% [u] = RotatingHillLShapeAdvec2D(u,FinalTime, G);
% 
% ua = Advec2Dsol(G.x,G.y,FinalTime,lam);
% err = ua - u; % compute point-wise error
% errL2 = zeros(G.K,1);
% for k = 1 : G.K
%     errL2(k) = err(:,k)'*diag(G.J(:,k))*G.MassMatrix*err(:,k);
% end
% errL2 = sqrt(sum(errL2)) % Global L^2-norm of error
% 
% figure
% for k = 1:G.K
%     T = delaunay(G.x(:,k),G.y(:,k));
%     trisurf(T,G.x(:,k),G.y(:,k),abs(u(:,k)))%-Advec2Dsol(x(:,k),y(:,k),FinalTime,lam))
%     hold on
% end


% Order of polymomials used for approximation 
NN = 3:4;
hh = 1/50;
E = zeros(length(NN), length(hh));
iterh = 0;
for hgrid = hh
    iterN = 0;
    iterh = iterh + 1;

    figure1 = figure;
    fd=@(p) drectangle(p,-1,1,-1,1);
    [p,G.EToV]=distmesh2d(fd,@huniform,hgrid,[-1,-1;1,1],[-1,-1;-1,1;1,-1;1,1]);
    close(figure1)

    G.VX = p(:,1)';
    G.VY = p(:,2)';
    
    G.EToV = Reorder(G.EToV, G.VX, G.VY);
    
    G.K = height(G.EToV);
    G.Nv = length(G.VX);
        
    fprintf("h = %d out of %d\n", iterh, length(hh));
    for N = NN
        G.N = N;
        iterN = iterN + 1;
        fprintf("\tN = %d out of %d\n", iterN, length(NN));
        
        % % Initialize solver and construct grid and metric
        G = StartUp2D(G);

        % Set initial conditions
        cx = 2*pi*G.y;
        cy = -2*pi*G.x;
        lam = 1/8;
        u = Advec2Dsol(G.x, G.y, 0, lam);

        % Solve Problem
        FinalTime = 1;
        [u] = RotatingHillLShapeAdvec2D(u,FinalTime, G);

        ua = Advec2Dsol(G.x,G.y,FinalTime,lam);
        err = ua - u; % compute point-wise error
        errL2 = zeros(G.K,1);
        for k = 1 : G.K
            errL2(k) = err(:,k)'*diag(G.J(:,k))*G.MassMatrix*err(:,k);
        end
        errL2 = sqrt(sum(errL2)); % Global L^2-norm of error
        E(iterN, iterh) = errL2;
        save('RotatingHillError.mat','E');
    end
end

%%
close all
figure
loglog(hh, E', '-o', LineWidth=2)
p = 1/2*round(2*(log(E(:,1))-log(E(:,end)))/(log(hh(1))-log(hh(end))));
xlabel('h')
legend(GetLegendLabels("N=",NN,"","  p=",p), location='best')
saveas(gca, "problem_5_RH_central_convergence.png")
