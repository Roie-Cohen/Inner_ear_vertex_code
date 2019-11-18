function area=cellarea(g,cnum)

area=0;
vidx = g.bonds(g.cells{cnum+1},1);
vert= getRelativePosition(g,vidx,cnum);
%%
vert = vert*g.scale; 
area = polyarea(vert(:,1),vert(:,2));
