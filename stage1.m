function g = stage1(g, ts, lat_ind)
% Gets an initial lattice 'g'.
% Runs stage 1 of the simulation for 'ts' time points.
% The function also saves the state of 'g' every few time points.
% The simulations run in parallel mode with the parameters:
% g.paras = [12 ; 0.3; 0.15; 0; 1; 1; 0.5];
% populations: 1 - general, 2 - SC, 3 - HC, 4 - Pillar, 5 - top border cell
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

savefold = 'stage1_timepoints'; % folder to save timepoints in

%% initial settings of the model and lattice

T1prob = 1; % probability for intercalation for small bonds
T1eps = 0.1; % bond length threshold
diff_reg = 2.8; % defines the region of differentiation above the PCs

g.transitionedBonds = [0, 0]; % keeps track of transitioned bonds. In the format of [bond_index, time_of_T1]
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

% weights: compressibility, tension, contractility, HCs repulsion, PCs pull, shear flow
g.paras = [12 ; 0.3; 0.15; 1; 1; 0.5]; 

g = LIdifferentiation(g, diff_reg);

timer = 0;
noise = rand(1, 2*length(g.verts));
g.globs = struct('timer',timer,'noise',noise);

%% run stage 1

save([savefold, '/lat(', num2str(lat_ind),')_step(0).mat'],'g');
for t =1:ts
    g = findTransitions(g,T1eps,T1prob,0.02);
    g = updateParameters(g);
    g = relaxLattice(g,5);
    
    g = killSmallCells(g,0.1);
    g = updateParameters(g);
    g = relaxLattice(g,5);
    
    g = LIdifferentiation(g, diff_reg);
    g = updateParameters(g);
    g = relaxLattice(g,10);
    
    if mod(t, 10) == 0
        save([savefold, '/lat(', num2str(lat_ind),')_step(',num2str(t/10),').mat'],'g');
    end
    
end


end