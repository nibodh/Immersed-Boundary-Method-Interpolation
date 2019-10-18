function [x, tf, side] = round_grid(m,dm,n,grid)

s = inputname(2);
if m < n
    side = 1;
    if strcmp(s,'dx')
        x = grid(1,ceil(m/dm)+1,1);
    elseif strcmp(s,'dy')
        x = grid(ceil(m/dm)+1,1,2);
    end
    if x < n
        tf = 1;
    else
        tf = 0;
    end
elseif m > n
    side = -1;
    if strcmp(s,'dx')
        x = grid(1,floor(m/dm)+1,1);
    elseif strcmp(s,'dy')
        x = grid(floor(m/dm)+1,1,2);
    end
    if x > n
        tf = 1;
    else
        tf = 0;
    end
else
    side = 0;
    tf = 0;
    if strcmp(s,'dx')
        x = grid(1,m/dm+1,1);
    elseif strcmp(s,'dy')
        x = grid(m/dm+1,1,2);
    end
end