clear all
close all

addpath 'Mesh2d v23'
addpath 'distmesh'
addpath '..\Codes1D'
addpath 'ServiceRoutines'
addpath '..\RiemannSolver'

clc

% Globals2D;
g=9.81;

% Polynomial order used for approximation 
S.N = 6;
S.hgrid=0.4;

sim='analyticCircle';
switch sim
    case {'circular'}
        [S.p,S.EToV,S.VX,S.VY,S.K,S.Nv,S.hgrid]=circularMesh2D(S.hgrid);
        BC  = @circularDamBC2D;
        IC = @circularDamIC2D;
        BC2  = @emptyDerivativeBC2D;
        FinalTime = 0.1;
        S=StartUp2D(S);
        S=circularBCtype(S);
        S.straight=1:S.K;
        S.curved=[];
    
    case {'wall'}
        [S.p,S.EToV,S.VX,S.VY,S.K,S.Nv,S.hgrid]=wallMesh2D(S.hgrid);
        BC  = @wallDamBC2D;
        IC = @wallDamIC2D;
        FinalTime = 7.2;
        S=StartUp2D(S);
        S=wallBCtype(S);
        BC2  = @wallDerivativeBC2D;
        S.straight=1:S.K;
        S.curved=[];

    case {'analytic'}
        [S.p,S.EToV,S.VX,S.VY,S.K,S.Nv,S.hgrid]=analyticMesh2D(S.hgrid);
        BC  = @analyticBC2D;
        IC = @analyticIC2D;
        BC2  = @emptyDerivativeBC2D;
        FinalTime = 0.01;
        S=StartUp2D(S);
        S=analyticBCtype(S);
        S=analyticSource2D(S);
        S.straight=1:S.K;
        S.curved=[];

    case {'analyticViscosity'}
        [S.p,S.EToV,S.VX,S.VY,S.K,S.Nv,S.hgrid]=analyticCircleMesh2D(S.hgrid);
        BC  = @analyticBC2D;
        IC = @analyticIC2D;
        BC2  = @analyticViscosityBC2D;
        FinalTime = 0.1;
        S=StartUp2D(S);
        S=analyticCircleBCtypeNoCurv(S);
%         S=analyticCircleBCtype(S);
        S=analyticSourceViscosity2D(S);
        S.straight=1:S.K;
        S.curved=[];


    case {'analyticCircle'}
        [S.p,S.EToV,S.VX,S.VY,S.K,S.Nv,S.hgrid]=analyticCircleMesh2D(S.hgrid);
        BC  = @analyticBC2D;
        IC = @analyticIC2D;
        BC2  = @emptyDerivativeBC2D;
        FinalTime = 0.1;
        S=StartUp2D(S);
        S=analyticCircleBCtypeNoCurv(S);
%         S=analyticCircleBCtype(S);
        S=analyticSource2D(S);
        S.straight=1:S.K;
        S.curved=[];

   case {'parabolic'}
        [S.p,S.EToV,S.VX,S.VY,S.K,S.Nv,S.hgrid]=parabolicMesh2D(S.hgrid);
        BC  = @parabolicBC2D;
        IC = @parabolicIC2D;
        BC2  = @parabolicDerivativeBC2D;
        FinalTime = 200;
        S=StartUp2D(S);
        S=parabolicBCtype(S);
        S=parabolicSource2D(S);
        S.straight=1:S.K;
        S.curved=[];

    case {'guus'}
        [S.p,S.EToV,S.VX,S.VY,S.K,S.Nv,S.hgrid]=guusMesh2D(S.hgrid);
        BC  = @guusBC2D;
        IC = @guusIC2D;
        BC2  = @guusDerivativeBC2D;
        FinalTime = 3;
        S=StartUp2D(S);
        S=guusBCtype(S);
end

figure
triplot(S.EToV,S.VX,S.VY)
hold on
plot(S.x(S.vmapO),S.y(S.vmapO),'*')


% axis equal
[TRI,xout,yout,interp,S] = plotStart(S.N,S);

figure
U=IC(S.x,S.y);
 PlotField2D(xout, yout, U(:,:,1),1,interp,TRI);
% axis equal

%% Solve Problem
figure(10)
U=IC(S.x,S.y);
fluxtype='HLLC';

type='CubatureLimit';
form='weak';
viscositytype='MDA';
tic
switch type
    case {'Cubature'}
        [U] = SWE2DCubature(U,FinalTime,BC,S,fluxtype);
    case {'CubatureLimit'}
        [U] = SWE2DCubatureLimit(U,FinalTime,BC,S,fluxtype);
    case {'Limit'}
        [U] = SWE2DLimit(U,FinalTime,BC,S,fluxtype,form);
    case {'CubatureViscosity'}
        [U,hsave] = SWE2DCubatureViscosity(U,FinalTime,BC,BC2,S,fluxtype,viscositytype);
    case {'ViscosityCubatureViscosity'}
        [U] = SWE2DViscosityCubatureViscosity(U,FinalTime,BC,BC2,S,fluxtype,viscositytype);
    case {'ViscosityTest'}
        [U] = SWE2DViscosityTest(U,FinalTime,BC,BC2,S,fluxtype,viscositytype);
end
toc
%%
% hf=cos(S.x-FinalTime).*cos(S.y-FinalTime)+2;
%      
% max(abs(hf-U(:,:,1)),[],'all')

%%
% FinalTime=10;
% [U] = SWE2DCubatureViscosity(U,FinalTime,BC,BC2,S,fluxtype);
% load('UsaveLimit.mat');
[TRI,xout,yout,interp,S] = plotStart(S.N,S);

figure
switch sim
    case {'circular'}
        PlotField2D(xout, yout, U(:,:,1),1,interp,TRI);
        xlim([-20 20])
        ylim([-20 20])
        zlim([0 3])
    
    case {'wall'}
        PlotField2D(xout, yout, U(:,:,1),1,interp,TRI);
        xlim([0 200])
        ylim([0 200])
        zlim([0 12])
        plotWall();
        view(40,35)

%         PlotAVG(xout,yout,U(:,:,1),S,1,interp,TRI)
%         xlim([0 200])
%         ylim([0 200])
%         zlim([0 12])
%         plotWall();
%         view(40,35)

    case {'analytic'}
        PlotField2D(xout, yout, U(:,:,1),1,interp,TRI);
        xlim([-2 2])
        ylim([-2 2])
        zlim([0 5])
    
     case {'analyticViscosity'}
        PlotField2D(xout, yout, U(:,:,1),1,interp,TRI);
        xlim([-2 2])
        ylim([-2 2])
        zlim([0 5])


    case {'analyticCircle'}
        PlotField2D(xout, yout, U(:,:,1),1,interp,TRI);
        xlim([-1 1])
        ylim([-1 1])
        zlim([0 5])

    case {'parabolic'}
        PlotField2D(xout, yout, U(:,:,1)+(8<= S.x & S.x<=12).*(0.2-0.05*(S.x-10).^2),1,interp,TRI);
        hold on
        PlotField2D(xout,yout, (8<= S.x & S.x<=12).*(0.2-0.05*(S.x-10).^2),0,interp,TRI);
        xlim([0 25])
        ylim([0 1])
        zlim([0 0.5])

    case {'guus'}
        PlotField2D(xout, yout, U(:,:,1),1,interp,TRI);
        hold on
        [X,Y,Z]=cylinder(1,100);
        surf(X,Y,3*Z,'FaceColor',[0.9290 0.6940 0.1250],'Edgecolor','None');
        xlim([-5 15])
        ylim([-10 10])
        zlim([0 4])

        figure
         PlotField2D(xout, yout, sqrt(U(:,:,1).^2+U(:,:,2).^2)./U(:,:,1)./sqrt(g*U(:,:,1)),1,interp,TRI);
        hold on
        [X,Y,Z]=cylinder(1,100);
        surf(X,Y,3*Z,'FaceColor',[0.9290 0.6940 0.1250],'Edgecolor','None');
        xlim([-5 15])
        ylim([-10 10])
        zlim([0 4])
end


set(gcf,'color','w');

xlabel('$x$','Interpreter', 'Latex', 'FontSize', 12)
ylabel('$y$','Interpreter', 'Latex', 'FontSize', 12)
zlabel('$h$','Interpreter', 'Latex', 'FontSize', 12)

% ax = gca;
% exportgraphics(ax,'wallMDA.pdf','Resolution',300)
% save('UsaveMDA.mat','U','S');

% exportgraphics(ax,'wallEV.pdf','Resolution',300)
% save('UsaveEV.mat','U');

% exportgraphics(ax,'wallDB.pdf','Resolution',300)
% save('UsaveDB.mat','U');

% exportgraphics(ax,'wallCubature.pdf','Resolution',300)
% save('UsaveCubature.mat','U');
% 
% exportgraphics(ax,'wallLimit.pdf','Resolution',300)
% save('UsaveLimit.mat','U');

%%
% Utemp=U(:,:,1)-(xout.^5+yout.^3+xout.^2.*yout.^3+xout.^2);
% Utemp(:,[1 3])=Utemp(:,[1 3])*1.1;
% Utemp=Utemp*0.5;
% figure
% PlotField2D(xout, yout, Utemp,1,interp,TRI([1:36+36*3],:));
% hold on
% triplot(S.EToV([1 2 3 4],:),S.VX,S.VY,'color',[0 0 0])
% xlabel('x')
% 
% idx=[1 7 28];
% for i=1:4
%     for j=1:3
%         plot3([S.VX(S.EToV(i,j)),S.VX(S.EToV(i,j))],[S.VY(S.EToV(i,j)),S.VY(S.EToV(i,j))],[0, Utemp(idx(j),i)],'Color',[0 0 0])
%     end
% end
% 
% box off
% axis off
% set(gcf,'color','w');
% view(-35,35)
% 
% ax = gca;
% exportgraphics(ax,'DGINTRO.pdf','Resolution',300)

%%
% load('guussave.mat')
% [TRI,xout,yout,interp,S] = plotStart(S.N,S);
% 
% figure
% PlotField2D(xout, yout, U(:,:,1),1,interp,TRI);
%         hold on
%         [X,Y,Z]=cylinder(1,100);
%         surf(X,Y,3*Z,'FaceColor',[0.9290 0.6940 0.1250],'Edgecolor','None');
%         xlim([-5 15])
%         ylim([-10 10])
%         zlim([0 4])

%%
% load('wallsave.mat')
% [TRI,xout,yout,interp,S] = plotStart(S.N,S);
% 
% figure
% PlotField2D(xout, yout, U(:,:,1),1,interp,TRI);
%         xlim([0 200])
%         ylim([0 200])
%         zlim([0 12])
%         plotWall();




%%
% ax = gca;
% exportgraphics(ax,'Froudmap.pdf','Resolution',300)


% save('hsaveEV.mat','hsave');

% load('hsaveh2.mat')

% %%
% figure
% PlotField2D(S.N, S.x, S.y, hsave(:,:,10),S,1);
% set(gcf,'renderer','zbuffer')
% xlim([0 200])
% ylim([0 200])
% zlim([0 12])
% plotWall();
%%


%%
% figure
% 
% xlim([0 200])
% ylim([0 200])
% zlim([0 12])
% plotWall();


%%
% figure
% axis padded
% for n = 1:1:size(hsave,3)-1
%       clf;
%       PlotAVG(xout,yout,hsave(:,:,n),S,1,interp,TRI)
%       hold on
%       xlim([0 200])
%       ylim([0 200])
%       zlim([0 12])
%       plotWall();
%       view(55,35)
%       set(gca,'visible','off')
%       set(gca,'xtick',[],'ytick',[])
%       set(gca,'xtick',[],'ytick',[])
%       set(gcf,'color','w');
% 
%       hold off
%       F = getframe();
%       im = frame2im(F); 
%       [imind,cm] = rgb2ind(im,256); 
%   if n == 1 
%       imwrite(imind,cm,'test5.gif','gif', 'Loopcount',inf,"DelayTime",0.1); 
%   else 
%       imwrite(imind,cm,'test5.gif','gif','WriteMode','append',"DelayTime",0.1); 
%   end 
% end


% figure
% axis padded
% for n = 1:6:size(hsave,3)-1
%       clf;
%       PlotField2D(S.N, S.x, S.y, hsave(:,:,n),S,1);
%       hold on
%       xlim([-5 15])
%       ylim([-10 10])
%       zlim([0 4])
%       [X,Y,Z]=cylinder(1,100);
%       surf(X,Y,3*Z,'FaceColor',[0.9290 0.6940 0.1250],'Edgecolor','None');
%       view(-40,40)
%       set(gca,'visible','off')
%       set(gca,'xtick',[],'ytick',[])
%       set(gca,'xtick',[],'ytick',[])
%       set(gcf,'color','w');
% 
%       hold off
% %       exportgraphics(gcf,'testAnimated12.gif','Append',true);
%       F = getframe();
%       im = frame2im(F); 
%       [imind,cm] = rgb2ind(im,256); 
%   if n == 1 
%       imwrite(imind,cm,'test.gif','gif', 'Loopcount',inf,"DelayTime",0.1); 
%   else 
%       imwrite(imind,cm,'test.gif','gif','WriteMode','append',"DelayTime",0.1); 
%   end 
% end
% 
% 
% %%
% % solNorm(U,S);
% 
% % figure
% % PlotField2D(xout, yout, U(:,:,1),1,interp,TRI);
% % %%
% 
% figure
% 
% 
% for i=0:0
%     load(append('test',num2str(i),'.mat'));
%     x=points(:,1);
%     y=points(:,2);
% 
%     PlotField2D(x, y, Density,1,1,tri+1)
%     hold on
% %     t=1.0;
% %     hf=cos(x-t).*cos(y-t)+2;
% % 
% %     max(abs(hf(:)-Density(:)))   
% %     figure
% %         PlotField2D(x, y, abs(hf'-Density),1,1,tri+1)
% end
%         hold on
%         [X,Y,Z]=cylinder(1,100);
%         surf(X,Y,5*Z,'FaceColor',[0.9290 0.6940 0.1250],'Edgecolor','None');
%         xlim([-5 25])
%         ylim([-10 10])
%         zlim([0 5])
% 
% % xlim([0 200])
% % ylim([0 200])
% % zlim([0 12])
% % plotWall();
% %%
% 
% figure
% 
% for i=0:0
%     load(append('test',num2str(i),'.mat'));
%     x=points(:,1);
%     y=points(:,2);
% 
% %     PlotField2D(x, y, sqrt(Velocity(:,1).^2+Velocity(:,2).^2)./sqrt(g*Density)',1,1,tri+1)
%     Fr=sqrt((Velocity(:,1)./Density').^2+(Velocity(:,2)./Density').^2)./sqrt(g*Density)';
% %     Fr(Fr>3)=3;
%     trisurf(tri+1, x(:), y(:), Fr,'Edgecolor','None')
%     colormap jet
%     colorbar
%     hold on
% 
% %     t=1.0;
% %     hf=cos(x-t).*cos(y-t)+2;
% % 
% %     max(abs(hf(:)-Density(:)))   
% %     figure
% %         PlotField2D(x, y, abs(hf'-Density),1,1,tri+1)
% end
%         hold on
%         [X,Y,Z]=cylinder(1,100);
%         surf(X,Y,3*Z,'FaceColor',[0.9290 0.6940 0.1250],'Edgecolor','None');
%         xlim([-3 7])
%         ylim([-3 3])
%         zlim([0 100])
%         view(2)
%         daspect([1 1 1])

%%
% figure
% plot(x,y,'*')
%%
% h=trisurf(double(tri)+1,x,y,Density)


% %%
% [TRI,xout,yout,interp,S] = plotStart(S.N+3,S);
% size(xout(:))
% h=interp*U(:,:,1);
% max(abs(sort(h(:))'-sort(double(Density))))
% 
% 
% max(abs(sort(xout(:))-sort(x)))
% max(abs(sort(yout(:))-sort(y)))

