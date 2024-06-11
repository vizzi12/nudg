clear all
close all
clc

figure
hold on
load('parabolic200.mat');
x=points(:,1);
y=points(:,2);
[x,ii]=sort(x(y< 1e-12));
z=(8<= x & x<=12).*(0.2-0.05*(x-10).^2);
h=Density(y< 1e-12)';

% load('ParabolicsaveMDA.mat');
% x=S.x;
% y=S.y;
% [x,ii]=sort(x(y< 1e-12));
% z=(8<= x & x<=12).*(0.2-0.05*(x-10).^2);
% h=U(:,:,1);
% h=h(y< 1e-12);


plot(x,h(ii)+z,'LineWidth',2,'color',[0 0.4470 0.7410 0.8]);
plot(x,z,'LineWidth',2,'color',[0.6250 0.3203 0.1758]);


x=0:0.01:8;
fx=0*x+0.4137357309;
plot(x,fx,'--','color',[0 0 0 0.8],'LineWidth',2);

x=8:0.01:10;
fx=(427716.*1i .* (x - 10) .* sqrt(0.30e2) .* sign((x - 10)) .* sqrt((18 .* 5 .^ (0.1e1 ./ 0.3e1) .* (x - 10) .^ 2 .* 327 .^ (0.2e1 ./ 0.3e1) + 324 .* 5 .^ (0.2e1 ./ 0.3e1) .* 327 .^ (0.1e1 ./ 0.3e1) + 109 .* (x - 10) .^ 4)) + (213858 .* 5 .^ (0.1e1 ./ 0.3e1) .* (x - 10) .^ 4 .* 327 .^ (0.2e1 ./ 0.3e1)) + (3849444 .* 5 .^ (0.2e1 ./ 0.3e1) .* (x - 10) .^ 2 .* 327 .^ (0.1e1 ./ 0.3e1)) + (1295029 .* x .^ 6) - (77701740 .* x .^ 5) + (1942543500 .* x .^ 4) - (25900580000 .* x .^ 3) + (194254350000 .* x .^ 2) - (777017400000 .* x) + 1294913516680) .^ (-0.1e1 ./ 0.3e1) .* ((6 .* 5 .^ (0.1e1 ./ 0.3e1) .* 327 .^ (0.2e1 ./ 0.3e1) + 109 .* (x - 10) .^ 2) .* (427716.*1i .* (x - 10) .* sqrt(0.30e2) .* sign((x - 10)) .* sqrt((18 .* 5 .^ (0.1e1 ./ 0.3e1) .* (x - 10) .^ 2 .* 327 .^ (0.2e1 ./ 0.3e1) + 324 .* 5 .^ (0.2e1 ./ 0.3e1) .* 327 .^ (0.1e1 ./ 0.3e1) + 109 .* (x - 10) .^ 4)) + (213858 .* 5 .^ (0.1e1 ./ 0.3e1) .* (x - 10) .^ 4 .* 327 .^ (0.2e1 ./ 0.3e1)) + (3849444 .* 5 .^ (0.2e1 ./ 0.3e1) .* (x - 10) .^ 2 .* 327 .^ (0.1e1 ./ 0.3e1)) + (1295029 .* x .^ 6) - (77701740 .* x .^ 5) + (1942543500 .* x .^ 4) - (25900580000 .* x .^ 3) + (194254350000 .* x .^ 2) - (777017400000 .* x) + 1294913516680) .^ (0.1e1 ./ 0.3e1) + (1308 .* 5 .^ (0.1e1 ./ 0.3e1) .* (x - 10) .^ 2 .* 327 .^ (0.2e1 ./ 0.3e1)) + (11881 .* x .^ 4) - (475240 .* x .^ 3) + (427716.*1i .* (x - 10) .* sqrt(0.30e2) .* sign((x - 10)) .* sqrt((18 .* 5 .^ (0.1e1 ./ 0.3e1) .* (x - 10) .^ 2 .* 327 .^ (0.2e1 ./ 0.3e1) + 324 .* 5 .^ (0.2e1 ./ 0.3e1) .* 327 .^ (0.1e1 ./ 0.3e1) + 109 .* (x - 10) .^ 4)) + (213858 .* 5 .^ (0.1e1 ./ 0.3e1) .* (x - 10) .^ 4 .* 327 .^ (0.2e1 ./ 0.3e1)) + (3849444 .* 5 .^ (0.2e1 ./ 0.3e1) .* (x - 10) .^ 2 .* 327 .^ (0.1e1 ./ 0.3e1)) + (1295029 .* x .^ 6) - (77701740 .* x .^ 5) + (1942543500 .* x .^ 4) - (25900580000 .* x .^ 3) + (194254350000 .* x .^ 2) - (777017400000 .* x) + 1294913516680) .^ (0.2e1 ./ 0.3e1) + (7128600 .* x .^ 2) + (11772 .* 5 .^ (0.2e1 ./ 0.3e1) .* 327 .^ (0.1e1 ./ 0.3e1)) - (47524000 .* x) + 118810000) ./ 6540;

z=2/10 - 1/20*(x - 10).^2;

plot(x,real(fx)+z,'--','color',[0 0 0 0.8],'LineWidth',2);

x=10:0.000001:11.6656;
fx=(0.109e3 ./ 0.120e3) .* (427716.*1i .* (x - 10) .* sign((x - 10)) .* sqrt(0.10e2) .* sqrt(0.3e1) .* sqrt((18 .* 5 .^ (0.1e1 ./ 0.3e1) .* (x - 10) .^ 2 .* 327 .^ (0.2e1 ./ 0.3e1) + 324 .* 5 .^ (0.2e1 ./ 0.3e1) .* 327 .^ (0.1e1 ./ 0.3e1) + 109 .* (x - 10) .^ 4)) + (213858 .* 5 .^ (0.1e1 ./ 0.3e1) .* (x - 10) .^ 4 .* 327 .^ (0.2e1 ./ 0.3e1)) + (3849444 .* 5 .^ (0.2e1 ./ 0.3e1) .* (x - 10) .^ 2 .* 327 .^ (0.1e1 ./ 0.3e1)) + (1295029 .* x .^ 6) - (77701740 .* x .^ 5) + (1942543500 .* x .^ 4) - (25900580000 .* x .^ 3) + (194254350000 .* x .^ 2) - (777017400000 .* x) + 1294913516680) .^ (-0.1e1 ./ 0.3e1) .* ((-0.1e1 ./ 0.11881e5.*1i .* sqrt(0.3e1) - (0.1e1 ./ 0.11881e5)) .* (427716.*1i .* (x - 10) .* sign((x - 10)) .* sqrt(0.10e2) .* sqrt(0.3e1) .* sqrt((18 .* 5 .^ (0.1e1 ./ 0.3e1) .* (x - 10) .^ 2 .* 327 .^ (0.2e1 ./ 0.3e1) + 324 .* 5 .^ (0.2e1 ./ 0.3e1) .* 327 .^ (0.1e1 ./ 0.3e1) + 109 .* (x - 10) .^ 4)) + (213858 .* 5 .^ (0.1e1 ./ 0.3e1) .* (x - 10) .^ 4 .* 327 .^ (0.2e1 ./ 0.3e1)) + (3849444 .* 5 .^ (0.2e1 ./ 0.3e1) .* (x - 10) .^ 2 .* 327 .^ (0.1e1 ./ 0.3e1)) + (1295029 .* x .^ 6) - (77701740 .* x .^ 5) + (1942543500 .* x .^ 4) - (25900580000 .* x .^ 3) + (194254350000 .* x .^ 2) - (777017400000 .* x) + 1294913516680) .^ (0.2e1 ./ 0.3e1) + (0.12e2 ./ 0.11881e5 .* (5 .^ (0.1e1 ./ 0.3e1)) .* (327 .^ (0.2e1 ./ 0.3e1)) + 0.2e1 ./ 0.109e3 .* ((x - 10) .^ 2)) .* (427716.*1i .* (x - 10) .* sign((x - 10)) .* sqrt(0.10e2) .* sqrt(0.3e1) .* sqrt((18 .* 5 .^ (0.1e1 ./ 0.3e1) .* (x - 10) .^ 2 .* 327 .^ (0.2e1 ./ 0.3e1) + 324 .* 5 .^ (0.2e1 ./ 0.3e1) .* 327 .^ (0.1e1 ./ 0.3e1) + 109 .* (x - 10) .^ 4)) + (213858 .* 5 .^ (0.1e1 ./ 0.3e1) .* (x - 10) .^ 4 .* 327 .^ (0.2e1 ./ 0.3e1)) + (3849444 .* 5 .^ (0.2e1 ./ 0.3e1) .* (x - 10) .^ 2 .* 327 .^ (0.1e1 ./ 0.3e1)) + (1295029 .* x .^ 6) - (77701740 .* x .^ 5) + (1942543500 .* x .^ 4) - (25900580000 .* x .^ 3) + (194254350000 .* x .^ 2) - (777017400000 .* x) + 1294913516680) .^ (0.1e1 ./ 0.3e1) + (0.36e2 ./ 0.109e3) .* ((x - 10) .^ 2) .* (1i .* (3 .^ (0.1e1 ./ 0.6e1)) .* (109 .^ (0.2e1 ./ 0.3e1)) - ((327 .^ (0.2e1 ./ 0.3e1)) ./ 0.3e1)) .* (5 .^ (0.1e1 ./ 0.3e1)) + 1i .* ((x - 10) .^ 4) .* sqrt(0.3e1) + (0.108e3 ./ 0.109e3.*1i .* (3 .^ (0.5e1 ./ 0.6e1)) .* (109 .^ (0.1e1 ./ 0.3e1)) - (0.108e3 ./ 0.109e3 .* (327 .^ (0.1e1 ./ 0.3e1)))) .* (5 .^ (0.2e1 ./ 0.3e1)) - ((x - 10) .^ 4));

z=2/10 - 1/20*(x - 10).^2;
% plot(x,real(fx)+z,'--','color',[0 0 0 0.8],'LineWidth',2);

x=11.6656:0.000001:12;

z2=2/10 - 1/20*(x - 10).^2;
fx2=x .^ 2 ./ 0.60e2 - x ./ 0.3e1 + 0.6785957e7 ./ 0.3956700e7 + (0.4348743025e10 .* x .^ 4 - 0.173949721000e12 .* x .^ 3 + 0.2634497078730e13 .* x .^ 2 - 0.17899997374600e14 .* x + (0.312436828833492232493e21 + 0.44522605742351046915e20 .* x .^ 2 - 0.182202918726222738300e21 .* x + 0.432664580055449775e18 .* x .^ 4 - 0.5835468850872991000e19 .* x .^ 3 + 0.286777858783625e15 .* x .^ 6 - 0.17206671527017500e17 .* x .^ 5 + 0.130571100e9 .* sqrt(-0.1720667152701750e16 .* x .^ 6 + 0.103240029162105000e18 .* x .^ 5 - 0.2595987480332698650e19 .* x .^ 4 + 0.35012813105237946000e20 .* x .^ 3 - 0.267135634454106281490e21 .* x .^ 2 + 0.1093217512357336429800e22 .* x - 0.1874774412310350284958e22)) .^ (0.2e1 ./ 0.3e1) + 0.46049212405849e14) .* (0.312436828833492232493e21 + 0.44522605742351046915e20 .* x .^ 2 - 0.182202918726222738300e21 .* x + 0.432664580055449775e18 .* x .^ 4 - 0.5835468850872991000e19 .* x .^ 3 + 0.286777858783625e15 .* x .^ 6 - 0.17206671527017500e17 .* x .^ 5 + 0.130571100e9 .* sqrt(-0.1720667152701750e16 .* x .^ 6 + 0.103240029162105000e18 .* x .^ 5 - 0.2595987480332698650e19 .* x .^ 4 + 0.35012813105237946000e20 .* x .^ 3 - 0.267135634454106281490e21 .* x .^ 2 + 0.1093217512357336429800e22 .* x - 0.1874774412310350284958e22)) .^ (-0.1e1 ./ 0.3e1) ./ 0.3956700e7;
plot([10:0.000001:11.6656,11.6656:0.000001:12],real([fx,fx2])+[z,z2],'--','color',[0 0 0 0.8],'LineWidth',2);





x=12:0.01:25;
fx=0*x+0.33;
plot(x,fx,'--','color',[0 0 0 0.8],'LineWidth',2)


xlim([0 25])
ylim([0 0.5])

xlabel('$x$','Interpreter', 'Latex', 'FontSize', 15)
ylabel('$h$','Interpreter', 'Latex', 'FontSize', 15)
legend('boxoff')  
box off
grid on
% axis square

legend('DG-FEM','','Analytic','Interpreter', 'Latex','FontSize', 15,'Location','northeast')

%%
% 
% ax=gca;
% exportgraphics(ax,'parabolicAnalyticZoom2.pdf','Resolution',300)






