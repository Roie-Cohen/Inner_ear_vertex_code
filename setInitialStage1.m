%% 
for i=31:100
    load(['random_lattices-211118\rand_lat(', num2str(i), ').mat'], 'g');
    g.populations(:) = 2;
%     for c=1:length(g.cells)-1
%         cy = mean(g.verts(g.bonds(g.cells{c+1},1),2));
%         if cy<-0.7 && cy >-1.3, g.populations(c) = 4; end
%     end
    g = selectCellType(g, 3, 4, 2);
    save(['initial_lats_stage1\s1ini_lat(', num2str(i), ').mat'], 'g');
    clf
end

%% flatten PCs
wb = waitbar(0);
for i=1:30
    load(['initial_lats_stage1\s1ini_lat(', num2str(i), ').mat'], 'g');
    
    for t=1:3
        g = updateTensions(g);
        g = relaxLattice(g, 10);
        g = findTransitions(g,0.1,1,0.02);
    end
   waitbar(i/30);
    save(['initial_lats_stage1\s1ini_latf(', num2str(i), ').mat'], 'g');
end
close(wb)