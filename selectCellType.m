function g = selectCellType(g, LCtype, RCtype, MCtype)
% selects types of cells by clicking on the lattice
% LCtype and RCtype is the type you choose by clicking right/left click
% esc to finish
figure(1)
LatticePresentation(g,0);

clist = find(~g.dead);
cpos = cellCenter(g, clist);
button = 1;
while button~=27
    [x,y,button] = ginput2(1);
    ds = cpos - [x y];
    Ds = sqrt(ds(:,1).^2 + ds(:,2).^2); % distances from every cell to the clicked point
    [m, mi] = min(Ds); % mi is the index in clist of the closest cell
    switch (button)
        case 1       % left mouse click
            g.populations(clist(mi)) = LCtype;
        case 2       % middle mouse click
            g.populations(clist(mi)) = MCtype;
        case 3       % right mouse click
            g.populations(clist(mi)) = RCtype;
    end
    LatticePresentation(g, 0);
end
close
end