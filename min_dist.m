function [x_min, y_min] = min_dist(xx,yy,dx,dy,grid,boundary,tolerance,x1,y1,x2,y2)
j = 0;    %%% counter for min distance
mark = zeros(length(xx));
dist = zeros(length(xx));
for i = 1:length(yy)
    mark(i) = cell_in(xx(i),yy(i),dx,dy,grid,boundary,tolerance);
    if mark(i) == 0
        j = j + 1;   %%% increasing counter for min distance on point found
        dist(i) = dist_ance([xx(i), yy(i)],[x1, y1]);
        if j == 1    %%% checking counter for min distance = 1
            distmin = dist(i);
            loc = i;
        elseif dist(i) < distmin
            distmin = dist(i);
            loc = i;
        end
    end
end
if j == 0
    y_min = y2;
    x_min = x2;
elseif j > 0
    y_min = yy(loc);
    x_min = xx(loc);
end
end