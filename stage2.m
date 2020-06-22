function g = stage2(g, ts, lat_ind)
% Gets a lattice 'g' at the final step of stage 1.
% Runs stage 2 of the simulation for 'ts' time points.
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

savefold = 'stage2_timepoints'; % folder to save timepoints in

%% initial settings of the model and lattice

T1prob = 1; % probability for intercalation for small bonds
T1eps = 0.1; % bond length threshold

g.transitionedBonds = [0, 0]; % keeps track of transitioned bonds. In the format of [bond_index, time_of_T1]
g.tension(2,2) = 3; % increasing the SC:SC tension

timer = 0;
noise = rand(1, 2*length(g.verts));
g.globs = struct('timer',timer,'noise',noise);

%% run stage 2
g = stage2differentiate(g);
save([savefold, '/lat(', num2str(lat_ind),')_step(0).mat'],'g');
for t =1:ts
    g = findTransitions(g,T1eps,T1prob,0.02);
    g = updateParameters(g);
    g = relaxLattice(g,5);
    
    g = killSmallCells(g,0.1);
    g = updateParameters(g);
    g = relaxLattice(g,5);
    
    if mod(t, 10) == 0
        save([savefold, '/lat(', num2str(lat_ind),')_step(',num2str(t/10),').mat'],'g');
    end
    
end


end