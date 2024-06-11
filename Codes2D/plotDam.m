clear all
close all
clc


load('wall72.mat');
x=points(:,1);
y=points(:,2);
idx=(100-1e-2<y & y< 100+1e-2);
[xsave,ii]=sort(x(idx));

h=Density(idx)';
hsave=h(ii);





hold on

load('WallSsave.mat')
load('Usavelimit.mat');
x=S.x;
y=S.y;
idx=(100-5e-1<y & y< 100+5e-1);
[x,ii]=sort(x(idx));
h=U(:,:,1);
h=h(idx);
h=h(ii);

plot(x,h,'-','LineWidth',2,'Color',[0 0.4470 0.7410 0.8])
plot(xsave,smooth(hsave),'--','LineWidth',2,'Color',[0 0 0 0.8])
xlim([140 200])
ylim([4.5 7])
xlabel('$x$', 'Interpreter', 'Latex', 'FontSize', 15)
ylabel('$h$', 'Interpreter', 'Latex','FontSize', 15)
legend('Limiter','Ref','Interpreter', 'Latex','FontSize', 15,'Location','northeast')
legend('boxoff')  
box off
axis square
grid on
set(gcf,'color','w');
ax=gca;
exportgraphics(ax,'limitZoom.pdf','Resolution',300)


figure
hold on

load('UsaveDB.mat');
h=U(:,:,1);
h=h(idx);
h=h(ii);

plot(x,h,'-','LineWidth',2,'Color',[0 0.4470 0.7410 0.8])
plot(xsave,smooth(hsave),'--','LineWidth',2,'Color',[0 0 0 0.8])
xlim([140 200])
ylim([4.5 7])
xlabel('$x$', 'Interpreter', 'Latex', 'FontSize', 15)
ylabel('$h$', 'Interpreter', 'Latex','FontSize', 15)
legend('DB','Ref','Interpreter', 'Latex','FontSize', 15,'Location','northeast')
legend('boxoff')  
box off
axis square
grid on
set(gcf,'color','w');
ax=gca;
exportgraphics(ax,'DBZoom.pdf','Resolution',300)


figure
hold on

load('UsaveEV.mat');
h=U(:,:,1);
h=h(idx);
h=h(ii);

plot(x,h,'-','LineWidth',2,'Color',[0 0.4470 0.7410 0.8])
plot(xsave,smooth(hsave),'--','LineWidth',2,'Color',[0 0 0 0.8])
xlim([140 200])
ylim([4.5 7])
xlabel('$x$', 'Interpreter', 'Latex', 'FontSize', 15)
ylabel('$h$', 'Interpreter', 'Latex','FontSize', 15)
legend('EV','Ref','Interpreter', 'Latex','FontSize', 15,'Location','northeast')
legend('boxoff')  
box off
axis square
grid on
set(gcf,'color','w');
ax=gca;
exportgraphics(ax,'EVZoom.pdf','Resolution',300)

figure
hold on

load('UsaveMDA.mat');
h=U(:,:,1);
h=h(idx);
h=h(ii);

plot(x,h,'-','LineWidth',2,'Color',[0 0.4470 0.7410 0.8])
plot(xsave,smooth(hsave),'--','LineWidth',2,'Color',[0 0 0 0.8])
xlim([140 200])
ylim([4.5 7])
xlabel('$x$', 'Interpreter', 'Latex', 'FontSize', 15)
ylabel('$h$', 'Interpreter', 'Latex','FontSize', 15)
legend('MDA','Ref','Interpreter', 'Latex','FontSize', 15,'Location','northeast')
legend('boxoff')  
box off
axis square
grid on
set(gcf,'color','w');
ax=gca;
exportgraphics(ax,'MDAZoom.pdf','Resolution',300)


figure
hold on

load('UsaveCubature.mat');
h=U(:,:,1);
h=h(idx);
h=h(ii);

plot(x,h,'-','LineWidth',2,'Color',[0 0.4470 0.7410 0.8])
plot(xsave,smooth(hsave),'--','LineWidth',2,'Color',[0 0 0 0.8])
xlim([140 200])
ylim([4.5 7])
xlabel('$x$', 'Interpreter', 'Latex', 'FontSize', 15)
ylabel('$h$', 'Interpreter', 'Latex','FontSize', 15)
legend('Naive','Ref','Interpreter', 'Latex','FontSize', 15,'Location','northeast')
legend('boxoff')  
box off
axis square
grid on
set(gcf,'color','w');
ax=gca;
exportgraphics(ax,'CubatureZoom.pdf','Resolution',300)
