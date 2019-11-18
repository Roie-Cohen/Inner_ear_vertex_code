function g = confineVerts(g)
% confines the vertices position to [-pi,pi]

nv = size(g.verts,1);
for i=1:nv
    for d=1:2
        if abs(g.verts(i,d)) > pi
            g.verts(i,d) = g.verts(i,d) - sign(g.verts(i,d))*2*pi;
        end
    end
end

end