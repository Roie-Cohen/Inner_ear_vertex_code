
warning off
global timer T1delay noise
wb = waitbar(0, 'Winter is coming...');
for lat_i=1:30
    
    load(['stage1_complete-211118\lat(',num2str(lat_i) ,')_step(20).mat'],'g');
    waitbar(lat_i/30);
    
    try
        T1prob = 1;
        T1eps = 0.1;
        
        T1delay = 0;
        timer = 0;
        noise = rand(1, 2*length(g.verts));
        
        for t =202:301
            
            g = findTransitions(g,T1eps,T1prob,0.02); %0.05, 0.02, 0.02
            g = updateParameters(g);
            g = relaxLattice(g,5); %50
            
            g = killSmallCells(g,0.1);
            g = updateParameters(g);
            g = relaxLattice(g,5); %200
            
            g = LIdifferentiation(g, 2.6);
            g = updateParameters(g);
            g = relaxLattice(g,10); %50
            
            if mod(t-1,10) == 0
                save(['stage1_complete-211118/lat(', num2str(lat_i),')_step(', num2str((t-1)/10) ,').mat'],'g');
            end
        end

    catch
        disp(['error in lattice: ', num2str(lat_i)]);
        continue;
    end
end

