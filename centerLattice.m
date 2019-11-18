function g = centerLattice(g)
    
    % calculating the central position of the lattice
    nc = length(g.cells)-1;
    centers = zeros(nc,2);
    for i=1:nc
        if ~g.dead(i)
            centers(i,:) = cellCenter(g,i);
        end
    end
    latt_center = mean(centers(~g.dead,:));
    
    % defining the lattice center as the origin
    nv = length(g.verts);
    g.verts = g.verts - [ones(nv,1)*latt_center(1) ones(nv,1)*latt_center(2) zeros(nv,1)];
    
end