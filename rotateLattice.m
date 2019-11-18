function g = rotateLattice(g, rad, isCentered)
% rotates the lattice around its center
    
    if ~isCentered
        g = centerLattice(g);
    end
    R = [cos(rad) -sin(rad); sin(rad) cos(rad)];
    g.verts(:,1:2) = g.verts(:,1:2)*R';

end