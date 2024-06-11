close all

addpath 'ServiceRoutines'
addpath '..\RiemannSolver'
% Driver script for solving the 1D advection equations
Globals1D;

g=9.81;

% % Order of polymomials used for approximation 

N = 8;
% % Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(0,50,50);
hgrid=VX(2)-VX(1);

% Initialize solver and construct grid and metric

StartUp1D;
% bx=zeros(size(x(1,:)));
% bx(x(2,:) > 10 & x(2,:) < 30)=-00.1;
% bx=ones(N+1,1).*bx;
bx=zeros(size(x));

h=3.5*ones(size(x(1,:)));
h(x(2,:) > 20)=1.25;
h=ones(N+1,1).*h;
u=zeros(size(x));
q=h.*u;

% h=0.33-(8<= x & x<=12).*(0.2-0.05*(x-10).^2);
% q=0.18*ones(size(x));
% % q=zeros(size(x));
% bx=(8<= x & x<=12).*(-0.10*x + 1.00);

size(x)
plot(x,h)
FinalTime = 2.5;
[h,q] = NSWE1D(h,q,FinalTime,bx);

x2=x;
h2=h;
%%
% figure
% plot(x(:),h(:)+(8<= x(:) & x(:)<=12).*(0.2-0.05*(x(:)-10).^2),'linewidth',2)


%%
figure
plot(x2(:),h2(:),'linewidth',2)


load dambreakdata
hold on
plot(x,h,'linewidth',2)
legend('DG-FEM','\texttt{dambreakdata}','location','northeast','Interpreter', 'Latex', 'FontSize', 15)
xlabel('$x$', 'Interpreter', 'Latex', 'FontSize', 15)
ylabel('$ h$', 'Interpreter', 'Latex','FontSize', 15)
ylim([0 4])

%
figure
plot(x2(:),q(:),'linewidth',2)
hold on
plot(x,hu,'linewidth',2)
legend('DG-FEM','\texttt{dambreakdata}','location','best','Interpreter', 'Latex', 'FontSize', 15)
xlabel('$x$', 'Interpreter', 'Latex', 'FontSize', 15)
ylabel('$ u$', 'Interpreter', 'Latex','FontSize', 15)
ylim([0 6])
