function runAlternativeModel(savefold, num_sim, tfinal)
% 22.05.2019. sim_date is a text variable includes the date of simulation
% runs the simulations for the alternative model in parallel mode

warning off

T1prob = 1;
T1eps = 0.1;
diff_reg = 3; % defines the region of differentiation above the PCs

gs = cell(num_sim, 1);
% set initial parameters for the cells
for lat_i=1:num_sim
    load(['initial_lats_stage1/s1ini_latf(', num2str(lat_i),').mat'],'g');
    g.transitionedBonds = [0, 0];
    g.area_feedback = 0; % no mechano-signaling
    g.fa = [1 1 1 1 1];  % General,SC,HCs,PC,top_border area ratio
    g = redistributeAreas(g);
    g.cmpr = ones(length(g.cells) - 1, 1); % compressibility of cells
    g.ctrc = ones(length(g.cells) - 1, 1); % contractility of cells
    g.centroid = cellCOM(g,1:length(g.cells)-1); % centroids of cells
    g.tension = [1  1  1  10 1 ; ... % tension of bond seperating populations n and m is g.tension(n,m)
                 1  3  1  10 10; ...
                 1  1  5  10 10; ...
                 10 10 10 1  1 ; ...
                 1  10 10 1  1 ];
    g.Gamma = [1 1 1 1 1]; % contractility of each cell type
    g.alpha = [1 1 1 1 1]; % compressibility of each cell type
    g.isAreaDiverg = 0; % 0 - no diverging term in the area energy. 1 - with diverging term.
    % weights: compressibility, tension, contractility, mechanical feedback, HCs repulsion, PCs pull, axial flow
    g.paras = [12 ; 0.3; 0.15; 0; 0; 0; 0]; % [12 ; 0.3; 0.15; 0; 1; 1; 0.5];
    
    g = LIdifferentiation(g, diff_reg);
    g = stage2differentiate(g);
    
    timer = 0;
    noise = rand(1, 2*length(g.verts));
    g.globs = struct('timer',timer,'noise',noise);
    
    gs{lat_i,1}=g;
    save([savefold, '/lat(', num2str(lat_i),')_step(0).mat'],'g');
end

for t =1:tfinal
    parfor lat_i=1:num_sim
        try
            
            g = gs{lat_i};
            
            g = findTransitions(g,T1eps,T1prob,0.02);
            g = updateParameters(g);
            g = relaxLattice(g,5);
            
            g = killSmallCells(g,0.1);
            g = updateParameters(g);
            g = relaxLattice(g,5);
            
            g = LIdifferentiation(g, diff_reg);
            g = updateParameters(g);
            g = relaxLattice(g,10);
            
            gs{lat_i}=g;
        catch
            disp(['error in lattice ', num2str(lat_i)]);
            gs{lat_i}=[];
            continue;
        end
    end
    
    % save the simulations in the current step
    for lat_i=1:num_sim
        if isempty(gs{lat_i}), continue; end
        g = gs{lat_i};
        save([savefold, '/lat(', num2str(lat_i),')_step(',num2str(t),').mat'],'g');
    end
    
    disp(['Total process: step number ', num2str(t), '. Percentage: ', num2str(round(100*t/tfinal, 2))]);
    
end

end


