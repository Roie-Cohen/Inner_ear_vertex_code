warning off
% parpool(8);
% global timer T1delay noise
% fileID = fopen('error_log.txt','a');
% savefold = '/a/home/cc/tree/taucc/students/physics/roiec/Sprinzak/vertex_code/stage1_complete_try';
savefold = 'stage1_repulsion';

ps = [0 0.2 0.5 1  1.5 3 5]; % parameters
for pi = 1:length(ps)
    
gs = cell(20,1);
for lat_i=1:20 
%     load(['/a/home/cc/tree/taucc/students/physics/roiec/Sprinzak/vertex_code/initial_lats_stage1/s1ini_latf(', num2str(lat_i),').mat'],'g');
    load(['stage1_complete-231218/lat(', num2str(lat_i),')_step(0).mat'],'g');
    gs{lat_i}=g;
end

parfor lat_i=1:20
    try
        
        g = gs{lat_i};
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
        % change in parameter:
        g.paras(5) = ps(pi);
        
        T1prob = 1;
        T1eps = 0.1;
        
        timer = 0;
        noise = rand(1, 2*length(g.verts));
        globs = struct('timer',timer,'noise',noise);
        for t =1:251 %1:251
            
            g = findTransitions(g,T1eps,T1prob,0.02,globs); %0.05, 0.02, 0.02
            g = updateParameters(g);
            [g, globs] = relaxLattice(g,5, globs);
            
            g = killSmallCells(g,0.1,globs);
            g = updateParameters(g);
            [g, globs] = relaxLattice(g,5, globs); 
            
            g = LIdifferentiation(g, 2.8);
            g = updateParameters(g);
            [g, globs] = relaxLattice(g,10, globs);
            
%             if mod(t-1,10) == 0
%                 save([savefold, '/lat(', num2str(lat_i),')_step(', num2str((t-1)/10) ,').mat'],'g');
%             end

        end
        gs{lat_i}=g;

    catch
%         fprintf(fileID,'\n lat=%1.2f t=%1.2f',[lat_i t]);
        gs{lat_i}=[];
        continue;
    end
end

for lat_i=1:20 
    if isempty(gs{lat_i}), continue; end
    g = gs{lat_i};
    save([savefold, '/lat(', num2str(lat_i),')_p-repulsion=',num2str(ps(pi)),'_final.mat'],'g');
end
% fclose(fileID);

end