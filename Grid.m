clear all
clc

%%%%%%%%%%%%%%%%%    grid   %%%%%%%%%%%%%%%%%%%%%         %%%% generic %%%%%
elementsx = 481;
elementsy = 289;
X = 20; %1.456832121;
Y = 12; %1.160353588;

dx = X/(elementsx-1);
dy = Y/(elementsy-1);
delta = (dx^2 + dy^2)^0.5;
[grid(:,:,1),grid(:,:,2)] = meshgrid(0:dx:X, 0:dy:Y);  %%%%%%   upside-down Y coordinate
grid(:,:,3) = 0;        %%%%%% 3rd dimensions: X coord, Y coord, epsilon
%plot(grid(:,:,1),grid(:,:,2),'o')

%%%%%%%%%%%%%%%%   geometry    %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   Geometry for Circle %%%%%%%%
R = 0.5;
h = 5; k = 6;
N = 160;

Geometry = zeros(N,6);
for i = 0:N
    Geometry(i+1,1) = h + R*cos(i*2*pi/N);% - pi/N);
    Geometry(i+1,2) = k - R*sin(i*2*pi/N);% - pi/N);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometry = textread('NACA651510.dat');
Geometry=unique(Geometry,'rows','stable');
N = length(Geometry);
plot(Geometry(:,1),Geometry(:,2),'o')
hold on;

%%%%%%%%%%%%%%%%%   identify boundary grid points    %%%%%%%%%%%%%%%%      %%%% generic %%%%%

for i = 1:N
    [flag,Geometry(i,3),Geometry(i,4),Geometry(i,5),Geometry(i,6),grid] =...
        incell(Geometry(i,1),Geometry(i,2),dx,dy,grid);
end

boundary_address = find(ismember(grid(:,:,3),0.5));        % finding positions of boundary points
boundary = grid(boundary_address(:));
boundary(:,2) = grid(boundary_address(:)+elementsx*elementsy);
plot(grid(boundary_address(:)),grid(boundary_address(:)+elementsx*elementsy), '+', 'Color','g')

%%%%%%%%%%%%%%%   identifying those points inside boundary   %%%%%%%%%%%%%
%%%%%%%%%%%%   generic    %%%%%%%55

in = inpoly(boundary,Geometry(:,1:2));
in_address = find(ismember(in,1));
plot(boundary(in_address(:),1),boundary(in_address(:),2),'x','Color','r')

%%%%%%%%%%%%%%   closest point %%%%%%%     %%%%% can be made faster  %%%%
% p = 1; % variable just to plot
for i = 1:length(in_address)
    for j = 1:N
        dist = ((Geometry(j,1)-boundary(in_address(i),1))^2 + (Geometry(j,2)-boundary(in_address(i),2))^2)^0.5;
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

for i = 1:length(in_address)
    k = boundary(in_address(i),3);
    for j = 3+1:3+k
    plot(Geometry(boundary(in_address(i),j),1),Geometry(boundary(in_address(i),j),2),'x','Color','r')
    end
end
%%%%%%%%%   finding reflection of point about boundary at closest points    %%%%
%%%%% generic %%%%%
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
        if m2 == 0
            if mod(y2,dy) == 0
                if x2 > x1
                    p1 = ceil(x1/dx);
                    if grid(1,p1+1,1) < x2
                        x2 = grid(1,p1+1,1);
                        overlap = 2;
                    end
                else
                    p1 = floor(x1/dx);
                    if grid(1,p1+1,1) > x2
                        x2 = grid(1,p1+1,1);
                        overlap = 2;
                    end
                end
            end
        elseif m == 0
            if mod(x2,dx) == 0
                if y2 > y1
                    q1 = ceil(y1/dy);
                    if grid(q1+1,1,2) < y2
                        y2 = grid(q1+1,1,2);
                        overlap = 2;
                    end
                else
                    q1 = floor(y1/dy);
                    if grid(q1+1,1,2) > y2
                        y2 = grid(q1+1,1,2);
                        overlap = 2;
                    end
                end
            end
        elseif m2 > 0            
            if x2 > x1                
                p1 = ceil(x1/dx);
                q1 = ceil(y1/dy);
                if and(grid(1,p1+1,1) < x2, grid(q1+1,1,2) < y2)                    
                    %check if lies on line
                    if ((grid(q1+1,1,2)-y1)/(grid(1,p1+1,1)-x1) == (y2-y1)/(x2-x1))                        
                        x2 = grid(1,p1+1,1);
                        y2 = grid(q1+1,1,2);
                        overlap = 2;
                    end
                end
            else
                p1 = floor(x1/dx);
                q1 = floor(y1/dy);
                if and(grid(1,p1+1,1) > x2, grid(q1+1,1,2) > y2)
                    %check if lies on line
                    if ((grid(q1+1,1,2)-y1)/(grid(1,p1+1,1)-x1) == (y2-y1)/(x2-x1))
                        x2 = grid(1,p1+1,1);
                        y2 = grid(q1+1,1,2);
                        overlap = 2;
                    end
                end
            end
        elseif m2 < 0
            if x2 < x1
                p1 = floor(x1/dx);
                q1 = ceil(y1/dy);
                if and(grid(1,p1+1,1) > x2, grid(q1+1,1,2) < y2)
                    %check if lies on line
                    if ((grid(q1+1,1,2)-y1)/(grid(1,p1+1,1)-x1) == (y2-y1)/(x2-x1))
                        x2 = grid(1,p1+1,1);
                        y2 = grid(q1+1,1,2);
                        overlap = 2;
                    end
                end                
            else
                p1 = ceil(x1/dx);
                q1 = floor(y1/dy);
                if and(grid(1,p1+1,1) < x2, grid(q1+1,1,2) > y2)
                    %check if lies on line
                    if ((grid(q1+1,1,2)-y1)/(grid(1,p1+1,1)-x1) == (y2-y1)/(x2-x1))
                        x2 = grid(1,p1+1,1);
                        y2 = grid(q1+1,1,2);
                        overlap = 2;
                    end
                end
            end
        end
    end
    points(i,:) = [x0 y0 x1 y1 x2 y2 overlap];
end
%%%%%% deleting overlapping points
short_points = find(ismember(points(:,7),2));
exact_points = find(ismember(points(:,7),1));
if isempty(exact_points) == 1
    %%%%%%%%%%%%%% insert checking for line overlap of band-boundary-forcing
    %%%%%%%%%%%%%% point and a grid point
    points(exact_points(:),3:4)=zeros(length(exact_points),2);
    l = 1;
    for i = 1:length(exact_points)
        [tf, found1] = ismemberf(points(exact_points(i),1:2),points(:,3:4),'row','tol',tolerance2);
        for j = 1:length(found1)
            if found1(j) <= length(points)
                found(l) = found1(j);
                l = l+1;
           end
        end
    end
    if l ~= 1
        found=unique(found);
        points(found(:),:)=[];
    end
end
plot(points(:,3), points(:,4), '+', 'Color', 'k');
plot(points(:,5), points(:,6), '*', 'Color', 'm');