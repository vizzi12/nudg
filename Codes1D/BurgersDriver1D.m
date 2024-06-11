% Driver script for solving the 1D Burgers equations
Globals1D;
close all
color=[	[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; [0.9290, 0.6940, 0.1250]];
% Order of polymomials used for approximation 

% Generate simple mesh
xL = -10; xR = 50;
N = 8;
[Nv, VX, K, EToV] = MeshGen1D(xL,xR,10);

% Initialize solver and construct grid and metric
StartUp1D;

for epsilon=[0.1 1 10]
    plot(x(:),1./(cosh(epsilon*(x(:)+5.0)).^2-0)+1.0,'LineWidth',2);
    hold on
end
xlabel('$x$', 'Interpreter', 'Latex', 'FontSize', 15)
legend('$\varepsilon=0.1$','$\varepsilon=1$','$\varepsilon=10$','$\mathcal{O}(h^{N+1})$','Interpreter', 'Latex','FontSize', 15)
ylabel('$u(x,0)$', 'Interpreter', 'Latex', 'FontSize', 15)
%%
% hplot=cell(3,1);
figure
% i=1;
% for N=[6, 10, 16]
    StartUp1D;
    % Set initial conditions
    epsilon = 1;
    u = 1./(cosh(epsilon*(x+5.0)-0).^2)+1.0;
    
    % Solve Problem
    FinalTime = 10;
    [u] = Burgers1D(u,epsilon,xL,xR,FinalTime);
    % [u] = Burgers1DSSP(u,epsilon,xL,xR,FinalTime);

%     ua = 1./(cosh(epsilon*(x+5.0)-FinalTime).^2)+1.0;
%     hold on
%     [u]  = Burgers1Dfilter(u,epsilon,xL,xR,FinalTime,16);
%     plot1=plot(x,u,'--','linewidth',2,'color',color(1,:));
% 
%     u = 1./(cosh(epsilon*(x+5.0)-0).^2)+1.0;
%     [u]  = Burgers1DSSP(u,epsilon,xL,xR,FinalTime);
%     plot2=plot(x,u,'--','linewidth',2,'color',color(2,:));

% end

%%
% legend([plot1(1),plot2(1)],'Interpolation','Exact','location','northwest','Interpreter', 'Latex', 'FontSize', 15)
% xlabel('$x$', 'Interpreter', 'Latex', 'FontSize', 15)
% ylabel('$u(x,50)$', 'Interpreter', 'Latex','FontSize', 15)
