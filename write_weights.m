function write_weights(g,fname)
n = length(g.cells)-1;
perim = zeros(n,1);
for i= 1:n;
  perim(i) = cellPerimeter(g,i);
end
use = find(g.bonds(:,1)~=0);
dif = g.verts(g.bonds(use,1),1:2)-g.verts(g.bonds(use,2),1:2);
dif = dif';

len = sum(dif .*dif)';
len = sqrt(len);
w1 = len ./perim(g.bonds(use,3));
w2 =  len ./perim(g.bonds(use,4));
res = [g.bonds(use,3) g.bonds(use,4) w1 w2];
dlmwrite(fname,res);

