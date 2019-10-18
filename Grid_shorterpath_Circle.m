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
[Grid(:,:,1),Grid(:,:,2)] = meshgrid(0:dx:X, 0:dy:Y);  %%%%%%   upside-down Y coordinate
Grid(:,:,3) = 0;        %%%%%% 3rd dimensions: X coord, Y coord, epsilon
%scatter(Grid(:,:,1),Grid(:,:,2),'o')

%%%%%%%%%%%%%%%   geometry    %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   Geometry for Circle %%%%%%%%
R = 0.5;
h = 5; k = 6;
N = 160;

Geometry = zeros(N,6);
for i = 0:N
    Geometry(i+1,1) = h + R*cos(i*2*pi/N);% - pi/N);
    Geometry(i+1,2) = k + R*sin(i*2*pi/N);% - pi/N);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometry = textread('circle.dat');
Geometry=unique(Geometry,'rows','stable');
N = length(Geometry);
figure(1)
plot(Geometry(:,1),Geometry(:,2),'LineWidth',1,'Color','k')
hold on;

%%%%%%%%%%%5    defining tolerance %%%%%%%%%
% vpa(dx)
% vpa(dy)
% vpa(((Geometry(1,1)-Geometry(2,1))^2 + (Geometry(1,2)-Geometry(2,2))^2)^0.5)
tolerance1 = 10e-11;
tolerance2 = 10e-11

%%%%%%%%%%%%%%%%%   identify boundary Grid points    %%%%%%%%%%%%%%%%      %%%% generic %%%%%

for i = 1:N
    [flag,Geometry(i,3),Geometry(i,4),Geometry(i,5),Geometry(i,6),Grid] =...
        incell(Geometry(i,1),Geometry(i,2),dx,dy,Grid);
end

boundary_address = find(ismember(Grid(:,:,3),0.5));        % finding positions of boundary points
boundary = Grid(boundary_address(:));
boundary(:,2) = Grid(boundary_address(:)+elementsx*elementsy);
scatter(Grid(boundary_address(:)),Grid(boundary_address(:)+elementsx*elementsy),16, 'x', 'MarkerEdgeColor','m','LineWidth',0.75)

%%%%%%%%%%%%%%%   identifying those points inside boundary   %%%%%%%%%%%%%
%%%%%%%%%%%%   generic    %%%%%%%55

in = inpoly(boundary,Geometry(:,1:2));
in_address = find(ismember(in,1));
scatter(boundary(in_address(:),1),boundary(in_address(:),2),16,'x','MarkerEdgeColor','b','LineWidth',0.75)

%%%%%%%%%%%%%%   closest point %%%%%%%
% p = 1; % variable just to scatter
for i = 1:length(in_address)
    for j = 1:N
        dist = dist_ance(Geometry(j,1:2),boundary(in_address(i),1:2));
        if j == 1
            distmin = dist;
            boundary(in_address(i),3) = 1;    % marking point having no equvidistant points
            boundary(in_address(i),4) = j;    % storing index in "Circle" of closest point
            continue
        end
        if dist < distmin
            distmin = dist;
            boundary(in_address(i),3) = 1;    % marking point having no equvidistant points
            boundary(in_address(i),4) = j;    % storing index in "Circle" of closest point
        elseif dist == distmin
            boundary(in_address(i),3) = boundary(in_address(i),3) + 1;
            boundary(in_address(i),3+boundary(in_address(i),3)) = j;    % storing index in "Circle" of closest point
        end
    end
end

%%%%%%%%%   finding reflection of point about boundary at closest points    %%%%
%%%%% generic %%%%%
overlap = 0;
points = zeros(length(in_address),7);
for i = 1:length(in_address)
    k = boundary(in_address(i),3);
    k1 = boundary(in_address(i),4);
    kp1 = k1 + 1;
    km1 = k1 - 1;
    if k1 == 1
        km1 = N;
        kp1 = k1 + 1;
    elseif k1 == N
        kp1 = 1;
        km1 = k1 - 1;
    end
    % first calc tangent of slope
    x0 = boundary(in_address(i),1);   % simplifying 
    y0 = boundary(in_address(i),2);
    overlap = 0;
    if k == 1
        xdiff = Geometry(kp1,1)-Geometry(km1,1);
        ydiff = Geometry(kp1,2)-Geometry(km1,2);
        px = Geometry(k1,1);
        py = Geometry(k1,2);
    elseif mod(k,2) == 0   %%% j-1 point is found
        px1 = Geometry(boundary(in_address(i),3+0.5*k),1);
        py1 = Geometry(boundary(in_address(i),3+0.5*k),2);
        px2 = Geometry(boundary(in_address(i),4+0.5*k),1);
        py2 = Geometry(boundary(in_address(i),4+0.5*k),2);
        xdiff = px1-px2;
        ydiff = py1-py2;
        px = (px1 + px2)/2;
        py = (py1 + py2)/2; 
    else
        xdiff = Geometry(boundary(in_address(i),3+floor(0.5*k)),1)...
            -Geometry(boundary(in_address(i),5+floor(0.5*k)),1);
        ydiff = Geometry(boundary(in_address(i),3+floor(0.5*k)),2)...
            -Geometry(boundary(in_address(i),5+floor(0.5*k)),2);     
        px = Geometry(boundary(in_address(i),4+floor(0.5*k)),1);
        py = Geometry(boundary(in_address(i),4+floor(0.5*k)),2);
    end
    
    % slope of line ------ determining surface and interpolant points
    if and(xdiff ~=0,ydiff ~= 0)
        m = ydiff/xdiff;
        m2 = -1/m;
        % calc which side
        syms x y
        eqn1 = m*x - y == m*px - py;
        eqn2 = m2*x - y == m2*x0 - y0;
        [y1, x1] = solve([eqn1, eqn2], y, x);
        if m*x0 - y0 + py - m*px > 0
            syms x y
            eqn1 = m*x - y == m*px - py - delta*(((m^2) +1)^0.5);
            [y2, x2] = solve([eqn1, eqn2], y, x);
        elseif m*x0 - y0 + py - m*px < 0
            syms x y
            eqn1 = m*x - y == m*px - py + delta*(((m^2) +1)^0.5);
            [y2, x2] = solve([eqn1, eqn2], y, x);
        else
            x2 = x0;
            y2 = y0;
            overlap = 1;  %overlap
        end        
    elseif xdiff ==0
        m2 = 0;
        y1 = y0;
        x1 = px;
        y2 = y0;
        if px > x0            
            x2 = px + delta;
        elseif px < x0
            x2 = px - delta;
        else
            x2 = x0;
            overlap = 1;   %overlap
        end
    elseif ydiff ==0
        m = 0;
        x1 = x0;
        y1 = py;
        x2 = x0;
        if py > y0
            y2 = py + delta;
        elseif py < y0
            y2 = py - delta;
        else
            y2 = y0;
            overlap = 1;    %overlap
        end
    end
    %%%% checking gird point along perpendicular
    if overlap == 0
        [x3, tf1, side1] = round_grid(x1,dx,x2,Grid);
        [y3, tf2, side2] = round_grid(y1,dy,y2,Grid);
        if m2 == 0
            if mod(y2,dy) == 0
               if tf1 == 1
                   x2 = x3;
                   overlap = 2;
               end
            end
        elseif m == 0
            if mod(x2,dx) == 0
                if tf2 == 1
                    y2 = y3;
                    overlap = 2;
                end
            end
        elseif and(tf1 == 1, tf2 == 1)        % 1 x and 1 y
            if abs((y3-y1)/(x3-x1)-(y2-y1)/(x2-x1)) < tolerance2   % checking slope point intersection
                x2 = x3;
                y2 = y3;
                overlap = 2;
            else
                if or(x3 + dx < x2, x3 + dx > x2)         % 2 x and 1 y3
                    x4 = x3 + side1*dx;
                    y_x3 = y1 + (x3-x1)*m2;
                    y_x4 = y1 + (x4-x1)*m2;
                    x_y3 = x1 + (y3-y1)/m2;
                    xx = [x3, x4, x_y3];
                    yy = [y_x3, y_x4, y3];
                elseif or(y3 + dy < y2, y3 + dy > y2)     % 2 y and 1 x3
                    y4 = y3 + side2*dy;
                    x_y3 = x1 + (y3-y1)/m2;
                    x_y4 = x1 + (y4-y1)/m2;
                    y_x3 = y1 + (x3-x1)*m2;
                    yy = [y3, y4, y_x3];
                    xx = [x_y3, x_y4, x3];
                else                    % 1 x3 and 1 y3
                    y_x3 = y1 + (x3-x1)*m2;
                    x_y3 = x1 + (y3-y1)/m2;
                    xx = [x3, x_y3];
                    yy = [y_x3, y3];
                end
                [x2, y2] = min_dist(xx,yy,dx,dy,Grid,boundary(in_address(:),1:2),tolerance2,x1,y1,x2,y2);
                %x2 = x2 + side1*tolerance2*cos(atan(m2));
                %y2 = y2 + side1*tolerance2*sin(atan(m2));
            end
        elseif or(tf1 == 1, tf2 == 1)
            if tf1 == 1
                if or(x3 + dx < x2, x3 + dx > x2)        % 2 x
                    x4 = x3 + side1*dx;                            
                    y_x3 = y1 + (x3-x1)*m2;
                    y_x4 = y1 + (x4-x1)*m2;
                    xx = [x3, x4];
                    yy = [y_x3, y_x4];
                else                    % 1 x
                    y_x3 = y1 + (x3-x1)*m2;
                    xx = x3;
                    yy = y_x3;
                end
            elseif tf2 == 1
                if or(y3 + dy < y2, y3 + dy > y2)         % 2 y
                    y4 = y3 + side2*dy;                            
                    x_y3 = x1 + (y3-y1)/m2;
                    x_y4 = x1 + (y4-y1)/m2;
                    yy = [y3, y4];
                    xx = [x_y3, x_y4];                        
                else                    % 1 y
                    x_y3 = x1 + (y3-y1)/m2;
                    xx = x_y3;
                    yy = y3;
                end
            end
            [x2, y2] = min_dist(xx,yy,dx,dy,Grid,boundary(in_address(:),1:2),tolerance2,x1,y1,x2,y2);
            %x2 = x2 + side1*tolerance2*cos(atan(m2));
            %y2 = y2 + side1*tolerance2*sin(atan(m2));
        end
    end
    points(i,:) = [overlap x0 y0 x1 y1 x2 y2];
end

%%%%%% deleting overlapping points
exact_points = find(ismember(points(:,1),1));
if isempty(exact_points) == 0
    %%%%%%%%%%%%%% insert checking for line overlap of band-boundary-forcing
    %%%%%%%%%%%%%% point and a Grid point
    points(exact_points(:),4:5)=zeros(length(exact_points),2);
    l = 1;
    for i = 1:length(exact_points)
        [tf, found1] = ismemberf(points(exact_points(i),2:3),points(:,4:5),'row','tol',tolerance2);
        points(exact_points(i),:)=[];        
        in_address(exact_points(:),:)=[];
    end
end
scatter(points(:,4), points(:,5),16, 's', 'MarkerFaceColor', 'c');
scatter(points(:,6), points(:,7),17, 'd', 'MarkerFaceColor', 'b');
for len = 1:length(points)
    plot([boundary(in_address(len),1),points(len,4),points(len,6)],[boundary(in_address(len),2),points(len,5),points(len,7)],'--','Color','k','LineWidth',0.25)
end
xlim([(4.5-2*dx) (5.5+2*dx)]);
ylim([(5.5-2*dy) (6.5+2*dy)]);
set(gca,'GridLineStyle','-','XColor', [0.9 0.9 0.9],'YColor', [0.9 0.9 0.9]...
    ,'XTick',0:dx:X,'YTick',0:dy:Y)
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.75 6]);
grid on;
h = legend('Boundary','Neighbor pts.','Forcing pts.','Surface pts.','Interpolation pts.','Location','best');
set(h,'box','off')
print -dpng Interpolation_Circle.png -r800;

fileID = fopen('XY_in_2.dat','w');
for ii = 1:length(points)
    fprintf(fileID,'%12.8f\t',points(ii,2:3));
    fprintf(fileID,'\r\n');
end
fclose(fileID);

fileID = fopen('XY_on_2.dat','w');
for ii = 1:length(points)
    fprintf(fileID,'%12.8f\t',points(ii,4:5));
    fprintf(fileID,'\r\n');
end
fclose(fileID);

fileID = fopen('XY_out_2.dat','w');
for ii = 1:length(points)
    fprintf(fileID,'%12.8f\t',points(ii,6:7));
    fprintf(fileID,'\r\n');
end
fclose(fileID);