function [x, tf] = round_grid(m,dm,n,grid)
if m < n
    x = ceil(m/dm);
    if x < n
        tf = 1;
    else
        tf = 0;
    end
elseif m > n
    x = floor(m/dm);
    if x > n
        tf = 1;
    else
        tf = 0;
    end
end