function Elength=cellPerimeter(g,i)

Elength=0;
vidx=g.bonds(g.cells{i+1},1);
vert = getRelativePosition(g,vidx,i);
nb=length(vidx);
for j=1:nb
 n = norm(vert(j,:)-vert(mod(j,nb)+1,:));
 Elength=Elength+n;
end
end