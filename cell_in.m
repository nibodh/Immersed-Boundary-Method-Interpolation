function mark = cell_in(x,y,dx,dy,grid,boundary,tolerance)
[flag,l1,m1,l2,m2] = incell(x,y,dx,dy,grid);
mark = 0;
if flag == 1
    l1 = grid(1,round(l1),1);
    m1 = grid(round(m1),1,2);
    if ismemberf([l1, m1],boundary,'rows','tol',tolerance) == 1
        mark = 1;
    end
elseif flag == 2
    l1 = grid(1,round(l1),1);
    m1 = grid(round(m1),1,2);
    l2 = grid(1,round(l2),1);
    if ismemberf([l1, m1],boundary,'rows','tol',tolerance) ||...
            ismemberf([l2, m1],boundary,'rows','tol',tolerance) == 1
        mark = 1;
    end
elseif flag == 3
    l1 = grid(1,round(l1),1);
    m1 = grid(round(m1),1,2);
    m2 = grid(round(m2),1,2);
    if ismemberf([l1, m1],boundary,'rows','tol',tolerance) ||...
            ismemberf([l1, m2],boundary,'rows','tol',tolerance) == 1
        mark = 1;
    end
elseif flag == 4
    l1 = grid(1,round(l1),1);
    m1 = grid(round(m1),1,2);
    l2 = grid(1,round(l2),1);
    m2 = grid(round(m2),1,2);
    if ismemberf([l1, m1],boundary,'rows','tol',tolerance) || ...
            ismemberf([l2, m1],boundary,'rows','tol',tolerance) || ...
            ismemberf([l1, m2],boundary,'rows','tol',tolerance)|| ...
            ismemberf([l2, m2],boundary,'rows','tol',tolerance) == 1
        mark = 1;
    end
end
end