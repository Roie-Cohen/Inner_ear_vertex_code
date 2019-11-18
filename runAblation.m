% using last frames in dtage 2 of 22.05.19

warning off
stage2_fold = 'stage2_2019-05-22'; % folder with stage1 simulations
s2tf = 100; % final step in stage1
% ablated_cells = [4 125 5 104 126 115 139 131 142 119 121 109 2 133];
% savefold = 'stage2_2019-05-22';
T1prob = 1;
T1eps = 0.1;
tfinal = 15; % simulation length


for lat_i = 1:50
    
    fname = [stage2_fold,'/lat(', num2str(lat_i),')_step(',num2str(s2tf),').mat'];
    if exist(fname, 'file')
        load([stage2_fold,'/lat(', num2str(lat_i),')_step(',num2str(s2tf),').mat'],'g');
    else
        disp([fname, ': not found']);
        continue;
    end
    
    nc = length(g.cells)-1;
    ab_cells_i = zeros(nc, 1);
    for c=1:nc
        if g.dead(c), continue; end
        npop = g.populations(g.bonds(g.cells{c+1},4)); % types of the neighbors of cell 'c'
        if sum(ismember([1 4 5 6], npop))==0 % if there are no neighbors other than HCs and SCs
            ab_cells_i(c) = 1;
        end
    end
    ablated_cells = find(ab_cells_i);
    
    for ai = 1:length(ablated_cells)

        close all
        load([stage2_fold,'/lat(', num2str(lat_i),')_step(',num2str(s2tf),').mat'],'g');
        timer = 0;
        
        noise = rand(1, 2*length(g.verts));
        g.globs = struct('timer',timer,'noise',noise);
        %  sim_folder = ['stage', num2str(stg),'_2019-05-22'];
        
        g.ablated_cell = ablated_cells(ai);
        
%          vid = VideoWriter(strcat('Simulations_vids/ablation-lat(',num2str(23),')-cell(',num2str(ablated_cells(ai)),').avi'));
%          vid.Quality = 100; % 100 for best
%          vid.FrameRate = 12;
%          open(vid);
        
        for t =1:tfinal
            
            g = findTransitions(g,T1eps,T1prob,0.02);
            g = updateParameters(g);
            g = relaxLattice(g,10);
            
            g = killSmallCells(g,0.15);
            g = updateParameters(g);
            g = relaxLattice(g,10);
            
            
        end
%         close(vid)
        save(['ablated_cells/ablation-lat(',num2str(lat_i),')-cell(',num2str(ablated_cells(ai)),').mat'],'g');
    end
    disp(['lattice ', num2str(lat_i), ' finished']);
end