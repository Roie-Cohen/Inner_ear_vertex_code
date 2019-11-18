function dE=dHarea(g,c)

A=cellarea(g,c);
A0=g.areas(c);
dE = zeros(2*length(g.verts),1);
vidx = g.bonds(g.cells{c+1},1);
vlist= getRelativePosition(g,vidx,c);
l = length(vidx);
for i=1:l,
    dE(2*vidx(i)-1) = vlist(mod(i,l)+1,2)-vlist(mod(i-2,l)+1,2);
    dE(2*vidx(i)) = vlist(mod(i-2,l)+1,1)-vlist(mod(i,l)+1,1);
end
% dE = -0.5*g.paras(1)*dE*A0*(A/A0 - 1 - A0/A);

if isfield(g,'populations')
    dE = -0.5*g.paras(1)*g.cmpr(c)*dE*A0*(A/A0 - 1 - g.isAreaDiverg*(A0/A)); 
else
    dE = -0.5*g.paras(1)*g.cmpr(c)*dE*A0*(A/A0 - 1);
end

end