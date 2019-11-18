function v=extractverts(g)
l = length(g.verts);
v=zeros(2*l,1);
v(1:2:2*l-1)=g.verts(:,1);
v(2:2:2*l)=g.verts(:,2);
% for i=1:length(g.verts)
%     v(2*i-1:2*i)=g.verts(i,1:2)';
% end