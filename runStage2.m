% 22.05.2019
% runs the simulations for stage2 in parallel mode with the parameters: 
% g.paras = [12 ; 0.3; 0.15; 0; 1; 1; 1]
% 

warning off
stage1_fold = 'stage1_2019-05-22'; % folder with stage1 simulations
s1tf = 400; % final step in stage1
savefold = 'stage2_2019-05-22';
num_sim = 50; % number of simulations
T1prob = 1;
T1eps = 0.1;
diff_reg = 2.8; % defines the region of differentiation above the PCs
tfinal = 100; % simulation length

gs = cell(num_sim, 1);
for lat_i=1:num_sim
    load([stage1_fold,'/lat(', num2str(lat_i),')_step(',num2str(s1tf),').mat'],'g');
    
    timer = 0;
    noise = rand(1, 2*length(g.verts));
    g.globs = struct('timer',timer,'noise',noise);
    g.tension(2,2) = 3; % increasing the SC:SC tension
    
        
    g = stage2differentiate(g);
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


