function runStage1(sim_date, num_sim, tfinal, param_name, param_values)
% 22.05.2019. sim_date is a text variable includes the date of simulation
% runs the simulations for stage1 in parallel mode with the parameters:
% g.paras = [12 ; 0.3; 0.15; 0; 1; 1; 0.5];
% populations: 1 - general, 2 - SC, 3 - HC, 4 - Pillar, 5 - top border cell
% if 'param_name' was given to the function, the simulation runs for
% different values of the parameter given in 'param_values'.
% 'gl-alpha' - the global weight of the area term
% 'gl-gamma' - the global weight of the tension term
% 'gl-Gama' - the global weight of the contractility term
% 'gl-eta' - the global weight of the shear term
% 'gl-zeta' - the global weight of the compression term
% 'gl-sigma' - the global weight of the repulsion term
% 'kappa' - the exponent of the repulsion term
% 'HC-Gamma' - the contractility of HCs
% 'HC-alpha' - the compressibility of HCs
% 'HCSCt' - the tension in HC:SC bonds

warning off
if nargin == 1
%     savefold = ['stage1_', sim_date];
    num_sim = 50; % number of simulations
    tfinal = 400; % simulation length
end

if nargin > 3
    savefold = ['stage1_', sim_date, '_', param_name];
    if ~exist(savefold, 'dir')
        mkdir(savefold);
    end
else
    param_values = 1;
    param_name = 'non';
    savefold = ['stage1_', sim_date];
end

T1prob = 1;
T1eps = 0.1;
diff_reg = 2.8; % defines the region of differentiation above the PCs

for pari = 1:length(param_values)
    
    gs = cell(num_sim, 1);
    add_t = ['_p=', num2str(param_values(pari))];
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
                    1  1  1  10 10; ...
                    1  1  5  10 10; ...
                    10 10 10 1  1 ; ...
                    1  10 10 1  1 ];
        g.Gamma = [1 1 2 1 1]; % contractility of each cell type
        g.alpha = [1 1 2 1 1]; % compressibility of each cell type
        g.isAreaDiverg = 0; % 0 - no diverging term in the area energy. 1 - with diverging term.
        % weights: compressibility, tension, contractility, mechanical feedback, HCs repulsion, PCs pull, axial flow
        g.paras = [12 ; 0.3; 0.15; 0; 1; 1; 0.5]; % [12 ; 0.3; 0.15; 0; 1; 1; 0.5];
      
        g = LIdifferentiation(g, diff_reg);
        
        timer = 0;
        noise = rand(1, 2*length(g.verts));
        g.globs = struct('timer',timer,'noise',noise);
        
        p = param_values(pari);
        switch param_name
            case 'gl-alpha'
                g.paras(1) = p;
            case 'gl-gamma'
                g.paras(2) = p;
            case 'gl-Gamma'
                g.paras(3) = p;
            case 'gl-eta'
                g.paras(7) = p;
            case 'gl-zeta'
                g.paras(6) = p;
            case 'gl-sigma'
                g.paras(5) = p;
            case 'kappa'
                g.paras(1) = p; % fix
            case 'HC-Gamma'
                g.Gamma(3) = p;
            case 'HC-alpha'
                g.alpha(3) = p;
            case 'HCSCt'
                g.tension(2,3) = p;
                g.tension(3,2) = p;
        end

        gs{lat_i,1}=g;
        save([savefold, '/lat(', num2str(lat_i),')',add_t,'_step(0).mat'],'g');
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
            save([savefold, '/lat(', num2str(lat_i),')',add_t,'_step(',num2str(t),').mat'],'g');
        end
        
        disp(['Total process: step number ', num2str(t), '. Percentage: ', num2str(round(100*t/tfinal, 2))]);
        
    end
    
end

end

