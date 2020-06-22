function g = relaxLattice(g,n)

for i=1:n
    g.globs.timer = g.globs.timer + 1;
    
    g = confineVerts(g);
    g = updateCellsCenter(g);
    
    ve = extractverts(g);
    dE=denergy(ve,g);
    noE = norm(dE);
    
    % add noise
    n_weight = 0.2*noE;
    dE(dE~=0) = dE(dE~=0) + g.globs.noise(dE~=0)*n_weight;
    if mod(g.globs.timer, 30)
        g.globs.noise = rand(1, 2*length(g.verts)) - 0.5;
    end
    
    if(noE>1E-10)
        
        dE = 0.05*dE/noE;
        ve = ve-dE';
        
        g_prev = g;
        g = insertverts(ve,g);
        g = ReposeLatticeByPCs(g, g_prev); % global reposition of the lattice
        
    end
    
end

end