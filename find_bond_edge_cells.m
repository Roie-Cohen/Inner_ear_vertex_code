function c = find_bond_edge_cells(g, bo)
% find the two cells on the egdes of the specific bond (not the two
% seperated by the bond)

bc = g.bonds(bo, 3); % one of the cells seperated by the bond
bordering_cell_bonds = g.cells{bc + 1}; % this cells' bonds
bo_ind = find(bordering_cell_bonds == bo); % the bond index in the cell bonds array
nb = length(bordering_cell_bonds);
next_bind = mod(bo_ind, nb) + 1;
prev_bind = mod(bo_ind - 2, nb) + 1;
next_cell = sum(g.bonds(bordering_cell_bonds(next_bind), 3:4)) - bc;
prev_cell = sum(g.bonds(bordering_cell_bonds(prev_bind), 3:4)) - bc;
c = [prev_cell next_cell];

end