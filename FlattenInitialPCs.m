warning off
global timer T1delay noise
fileID = fopen('error_log.txt','a');
wb = waitbar(0, 'Winter is coming...');
for lat_i=1:50
    
    load(['initial_lats_stage1\s1ini_lat(', num2str(lat_i),').mat'],'g');
    waitbar(lat_i/50);
    
    try
        g.transitionedBonds = [0, 0];
        g.area_feedback = 0;
        g.fa = [1 1 0.7 1.2 1 1];  % Hensen:SC:HCs:pillar area ratio
        g = redistributeAreas(g);
        g.cmpr = ones(length(g.cells) - 1, 1); % compressibility of cells
        g.ctrc = ones(length(g.cells) - 1, 1); % contractility of cells
        g.area_pert = zeros(length(g.cells) - 1, 1); % logial array, if a cell's area is perturbed
        g.centroid = zeros(length(g.cells) - 1, 2); % centroids of cells
%         g = LIdifferentiation(g, 2.6);
        g.paras = [6 ; 0.3; 0.15; 0; 1; 0.5; 0.4];
        T1prob = 1;
        T1eps = 0.1;
        
        T1delay = 0;
        timer = 0;
        noise = rand(1, 2*length(g.verts));
        
        for t =1:10%251
            
            g = findTransitions(g,T1eps,T1prob,0.02); %0.05, 0.02, 0.02
            g = updateParameters(g);
            g = relaxLattice(g,5);
            
            g = killSmallCells(g,0.1);
            g = updateParameters(g);
            g = relaxLattice(g,5); 
            
%             g = LIdifferentiation(g, 2.8);
%             g = updateParameters(g);
%             g = relaxLattice(g,10); %50
            
%             if mod(t-1,10) == 0
%                 save(['stage1_complete-131218/lat(', num2str(lat_i),')_step(', num2str((t-1)/10) ,').mat'],'g');
%             end
        save(['initial_lats_stage1\s1ini_latf(', num2str(lat_i),').mat'],'g');
        end

    catch
        fprintf(fileID,'\n lat=%1.2f t=%1.2f',[lat_i t]);
        continue;
    end
end
fclose(fileID);
close(wb)