%JK plot
clear all
clc

%%%%%%%%%%%%%%%%%    Grid   %%%%%%%%%%%%%%%%%%%%%         %%%% generic %%%%%
elementsx = 401
elementsy = 241
X = 20 %1.456832121;
Y = 12 %1.160353588;

dx = X/(elementsx-1);
dy = Y/(elementsy-1);
delta = (dx^2 + dy^2)^0.5;

Geometry = textread('circle.dat');
figure(1)
plot(Geometry(:,1),Geometry(:,2),'LineWidth',1,'Color','k')
hold on;

In = textread('in.dat');
scatter(In(:,1),In(:,2),16,'x','MarkerEdgeColor','b','LineWidth',0.75)

points = textread('inBSI.dat');
scatter(points(:,3), points(:,4),16, 's', 'MarkerFaceColor', 'c');
scatter(points(:,5), points(:,6),17, 'd', 'MarkerFaceColor', 'b');
xlim([(4.5-2*dx) (5.5+2*dx)]);
ylim([(5.5-2*dy) (6.5+2*dy)]);
set(gca,'GridLineStyle','-','XColor', [0.9 0.9 0.9],'YColor', [0.9 0.9 0.9]...
    ,'XTick',0:dx:X,'YTick',0:dy:Y)
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.75 6]);
grid on;
h = legend('Boundary','Neighbor pts.','Forcing pts.','Surface pts.','Interpolation pts.','Location','best');
set(h,'box','off')
print -dpng JK.png -r800;