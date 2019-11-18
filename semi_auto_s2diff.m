function g = semi_auto_s2diff(g)
% semi-automated top-border differentiation of a lattice
g = stage2differentiate(g);
g = selectCellType(g,2,5,1);
end