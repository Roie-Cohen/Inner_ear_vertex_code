function gout=insertverts(v,gin)

gout=gin;
gout.verts(:,1:2) = reshape(v,[2 length(gout.verts)])';
% for i=1:length(gout.verts)
%     gout.verts(i,1:2)=v(2*i-1:2*i)';
% end