function g=rescale(gin)
g = gin;
xmin = min(g.verts(:,1));
xmax = max(g.verts(:,1));
ymin = min(g.verts(:,2));
ymax = max(g.verts(:,2));
xmean = mean(g.verts(:,1));
ymean = mean(g.verts(:,2));

g.verts(:,1) = g.verts(:,1)-xmean;
g.verts(:,2) = g.verts(:,2)-ymean;
end