function [len center]= getBoundaryLength(g,idx,cel)
if(nargin==2)
   cel = g.bonds(idx,3); 
end
vidx1 = g.bonds(idx,1);
vidx2 = g.bonds(idx,2);
pos = getRelativePosition(g,[vidx1 vidx2]',cel);
pos1 = pos(1,:);
pos2 = pos(2,:);
d = (pos2 -pos1)*g.scale; 
d = d .*d;
len = sqrt(sum(d,2));
center = 0.5*(pos1+pos2);
end