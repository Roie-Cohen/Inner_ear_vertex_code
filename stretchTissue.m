function g = stretchTissue(g, dx)
% works only for 12x12 hexagonal non-periodic lattice

edge_x_pos = [9.958, 10.82]; % need to improve this part
edge_verts = g.bonds(g.cells{1} ,1);
st_verts = edge_verts( abs(g.verts(edge_verts,1)) > edge_x_pos(1) );
g.verts(st_verts, 1) = sign(g.verts(st_verts, 1)).*(edge_x_pos(2) + dx);
g.fixed_verts = [st_verts, g.verts(st_verts, 1), g.verts(st_verts, 2)];

end