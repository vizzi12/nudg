clear all
close all
clc
figure
hold on
c=[[0 0.4470 0.7410 0.5]; [0.8500 0.3250 0.0980 0.5]; [0.6350 0.0780 0.1840 0.5]; [0.4660 0.6740 0.1880 0.5]];
for i=0:3
    load(append('test',num2str(i),'.mat'));

    x=points(:,1);
    y=points(:,2);

    triplot(double(tri)+1,x,y,'LineWidth',2,'color',c(i+1,:))
end
axis off

ax = gca;
exportgraphics(ax,'meshHalo.pdf','Resolution',300)
