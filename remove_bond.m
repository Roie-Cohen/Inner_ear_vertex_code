function g = remove_bond(g,bond)
g.bonds(bond,:)=0;
for c = 1:length(g.cells)-1
    g.cells{c+1} = g.cells{c+1}(~ismember(g.cells{c+1},bond));
end
