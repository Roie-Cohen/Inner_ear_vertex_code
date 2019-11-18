function [area, rel] = estimatedArea(g, c)
% estimation of the real area of HCs and SCs considering curvature.
% R the radius of curvature.
% add - indicates if the add or reduce the curvature areas. add=+1 or -1
% also returns the relative area between the real area and polygonal area.
A0 = cellarea(g, c);
area = A0;
vidx = g.bonds(g.cells{c+1},1);
verts = getRelativePosition(g, vidx);
add = 1;
if g.populations(c) == 2, add = -1; end
for i=1:length(vidx)
    adj_cells = g.bonds(g.cells{c+1}(i),3:4);
    p = g.populations(adj_cells);
    if sum(ismember([2 3], p)) ~= 2, continue; end % only HC:SC bonds
    hci = adj_cells(p==3); % the index of the HC
    A_hc = cellarea(g,hci); % the area of the HC
    sides_hc = length(g.cells{hci+1});
    r_fac = 2/(sides_hc*sin(2*pi/sides_hc)); % see methods
    R = sqrt(r_fac*A_hc); 
    L = norm(verts(i,:) - verts(mod(i,length(vidx))+1,:)); % length of the bond
    h = sqrt(R^2-0.25*L^2); % distance from the circle center to the bond
    alpha = 2*asin(0.5*L/R); % angle of the bond in the circle
    segment = 0.5*alpha*R^2-0.5*h*L;
    area = area + add*segment;
end
rel = area/A0;
end