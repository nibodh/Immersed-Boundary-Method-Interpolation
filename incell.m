function [flag,x1,y1,x2,y2,grid] = incell(x,y,dx,dy,grid)

    xrem = mod(double(x),dx);
    yrem = mod(double(y),dy);
    if and(xrem == 0, yrem == 0)
        flag = 1;
        x1 = x/dx + 1;    %%%%%% 2nd dimension X coord,, Y coord, cell
        y1 = y/dy + 1;    %%%% + 1 because grid starts with 0 - not 1
        x2 = 0;
        y2 = 0;
        grid(y1, x1, 3) = 0.5;
    elseif or(xrem == 0, yrem == 0)
        if yrem == 0
            flag = 2;
            x1 = floor(x/dx) + 1;
            y1 = y/dy + 1;
            x2 = ceil(x/dx) + 1;
            y2 = 0;
            grid(round(y1), round(x1), 3) = 0.5;
            grid(round(y1), round(x2), 3) = 0.5;
        elseif xrem == 0
            flag = 3;
            x1 = x/dx + 1;
            y1 = floor(y/dy) + 1;
            x2 = 0;
            y2 = ceil(y/dy) + 1;
            grid(round(y1), round(x1), 3) = 0.5;
            grid(round(y2), round(x1), 3) = 0.5;
        end
    else
        flag = 4;
        x1 = floor(x/dx) + 1;
        y1 = floor(y/dy) + 1;
        x2 = ceil(x/dx) + 1;
        y2 = ceil(y/dy) + 1;
        grid(round(y1), round(x1), 3) = 0.5;
        grid(round(y1), round(x2), 3) = 0.5;
        grid(round(y2), round(x1), 3) = 0.5;
        grid(round(y2), round(x2), 3) = 0.5;
    end
end