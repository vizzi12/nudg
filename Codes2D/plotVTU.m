
addpath 'Mesh2d v23'
addpath 'distmesh'
addpath '..\Codes1D'
addpath 'ServiceRoutines'
addpath '..\RiemannSolver'

sim='guus';
figure
for i=0:0
    load(append('test',num2str(i),'.mat'));
%     load('wall72.mat');
%     load('circularDam47.mat');
%     load('parabolic200.mat');
%     load(append('circularTest',num2str(i),'.mat'));

    x=points(:,1);
    y=points(:,2);
    switch sim
      case {'basic'}
%         [C,ia,ic] = uniquetol([x,y],'byrows',true);
%         a_counts = accumarray(ic,1);
% 
%         Density = accumarray(ic,Density)./a_counts;
%         Density = Density(ic);
        PlotField2D(x, y, Density,1,1,tri+1);
        hold on
        t=0.1;

        hf=cos(x-t).*cos(y-t)+2;
        max(abs(hf(:)-Density(:)))
      case {'circularDam'}
        [C,ia,ic] = uniquetol([x,y],'byrows',true);
        a_counts = accumarray(ic,1);
        Density = accumarray(ic,Density)./a_counts;
        Density = Density(ic);
        PlotField2D(x, y, Density,1,1,tri+1);
        hold on
        xlim([-20 20])
        ylim([-20 20])
        view(45,80)
%         zlim([0 1])

      case {'guus'}
        PlotField2D(x, y, Density,1,1,tri+1)
        hold on
        [X,Y,Z]=cylinder(1,100);
        surf(X,Y,5*Z,'FaceColor',[0.9290 0.6940 0.1250],'Edgecolor','None');
        xlim([-5 15])
        ylim([-10 10])
        zlim([0 5])

        figure
        g=9.81;
        Fr=sqrt((Velocity(:,1)./Density').^2+(Velocity(:,2)./Density').^2)./sqrt(g*Density)';
%         Fr=sqrt((Velocity(:,1)./Density').^2+(Velocity(:,2)./Density').^2);%./sqrt(g*Density)';
        trisurf(tri+1, x(:), y(:), Fr,'Edgecolor','None')
        colormap jet
        colorbar
        xlim([-3 7])
        ylim([-5 5])
        zlim([0 100])

        view(2)
      case {'wall'}
 
        
        [C,ia,ic] = uniquetol([x,y],'byrows',true);
        a_counts = accumarray(ic,1);

        Density = accumarray(ic,Density)./a_counts;
        Density = Density(ic);

        PlotField2D(x, y, Density,1,1,tri+1)
        hold on
        xlim([0 200])
        ylim([0 200])
        zlim([0 12])
        plotWall();
        view(40,35)
        set(gcf,'color','w');
      case {'parabolic'}
        [C,ia,ic] = uniquetol([x,y],'byrows',true);
        a_counts = accumarray(ic,1);

        Density = accumarray(ic,Density)./a_counts;
        Density = Density(ic);
        PlotField2D(x,y, (8<= x & x<=12).*(0.2-0.05*(x-10).^2),0,1,tri+1);
        hold on
        PlotField2D(x, y, Density+(8<= x & x<=12).*(0.2-0.05*(x-10).^2),1,1,tri+1);

        xlim([0 25])
        ylim([0 1])
        zlim([0 0.5])
        view(-15,20)

        
    end
    xlabel('$x$','Interpreter', 'Latex', 'FontSize', 12)
    ylabel('$y$','Interpreter', 'Latex', 'FontSize', 12)
    zlabel('$h$','Interpreter', 'Latex', 'FontSize', 12)
end

set(gcf,'color','w');

%%
% ax=gca;
% exportgraphics(ax,'parabolic.pdf','Resolution',300)