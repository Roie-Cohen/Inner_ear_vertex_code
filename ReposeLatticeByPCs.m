function g = ReposeLatticeByPCs(g, g_prev)
% repositions the lattice according to previous position of the PCs

PCs = find(g.populations == 4);
if isempty(PCs), return; end
PCs_pos = cellCOM(g, PCs);
current_pos = mean(PCs_pos,1);

% previous PCs position
prev_PCs = find(g_prev.populations == 4);
prev_PCs_pos = cellCOM(g_prev, prev_PCs);
prev_pos = mean(prev_PCs_pos,1);

dpos = prev_pos - current_pos;

g.verts(:,1:2) = g.verts(:,1:2) + dpos;

end