function plotSol(h,S,varargin)
    Globals2D;
    if(~isempty(varargin))
        figure(varargin{1});
    else
        figure
    end
    
    for i = 1:S.K
        tri = delaunay(S.x(:,i),S.y(:,i));
        trisurf(tri,S.x(:,i),S.y(:,i),h(:,i),'FaceColor',[0.3010 0.7450 0.9330],'Edgecolor','None','Edgecolor','None')
        hold on
    end
    material shiny,    lighting gouraud 
    camlight headlight,
    hold off
end