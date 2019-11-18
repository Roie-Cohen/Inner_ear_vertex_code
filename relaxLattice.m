function g = relaxLattice(g,n,vid)

for i=1:n
    g.globs.timer = g.globs.timer + 1;
    
    g = confineVerts(g);
    if g.populations(1)~=0, g = updateCellsCenter(g); end
    ve = extractverts(g);
    dE=denergy(ve,g);
    noE = norm(dE);

    dE(dE~=0) = dE(dE~=0) + g.globs.noise(dE~=0)*0.2*noE;
    if mod(g.globs.timer, 30)
        g.globs.noise = rand(1, 2*length(g.verts));
    end

    if(noE>1E-10),
        
        dE = 0.05*dE/noE; 
        ve = ve-dE';
        g_prev = g;
        g = insertverts(ve,g);
        if g.populations(1)~=0, g = ReposeLatticeByPCs(g, g_prev); end
        
        if nargin == 3
            % capturing a frame every nv steps
            nv = 5;
            if (mod(i,nv)==0 || g.globs.timer==1)
                figure(6);
                LatticePresentation(g,0, 6);
                if g.bc == 1
                    lm = 3.5;
                else
                    lm = 10;
                end
                xlim([-lm lm]);
                ylim([-lm lm]);
                frame = getframe;
                writeVideo(vid,frame);
            end
        end
    end
    
end

end