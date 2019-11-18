
warning off
global ctypes timer T1delay noise
load('lattices/LIlat_autodiff1','g');
g.transitionedBonds = [0, 0];
g.area_feedback = 0;
g.fa = [1 1 0.7 1.2 1 1];  % Hensen:SC:HCs:pillar area ratio
g = redistributeAreas(g);
g.cmpr = ones(length(g.cells) - 1, 1); % compressibility of cells
g.ctrc = ones(length(g.cells) - 1, 1); % contractility of cells
g.area_pert = zeros(length(g.cells) - 1, 1); % logial array, if a cell's area is perturbed
g.centroid = zeros(length(g.cells) - 1, 2); % centroids of cells
g = LIdifferentiation(g, 2.4);
g1 = g;

T1prob = 1;
T1eps = 0.1;

fileID = fopen('error_log.txt','a');
wb = waitbar(0, 'Winter is coming...');
counter = 1;
for p1=6:3:12
    for p2=0.1:0.1:0.3
        for p3=0.15:0.05:0.25
            for p6=0.5:0.25:1
                for p7=0.2:0.2:0.4
                    try
                    g = g1;
                    ctypes = [];
                    T1delay = 0;
                    timer = 0;
                    % parameters: area, tension, perimeter, mechanical feedback, HCs vdv, PCs
                    % pull, axial flow
                    g.paras = [p1 ; p2; p3; 0; 1; p6; p7];
                    fname = ['simulations-191118/lat_step(0)_comp(', num2str(p1), ')_tens(', num2str(p2), ')_cont(', num2str(p3), ')_pp(', num2str(p6), ')_sh(', num2str(p7), ').mat'];
                    if isfile(fname), continue; end
                    noise = rand(1, 2*length(g.verts));
                    for t =1:201
                        
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
                            save(['simulations-191118/lat_step(', num2str((t-1)/10), ')_comp(', num2str(p1), ')_tens(', num2str(p2), ')_cont(', num2str(p3), ')_pp(', num2str(p6), ')_sh(', num2str(p7), ').mat'],'g');
                        end
                    end
                    waitbar(counter/162);
                    counter = counter + 1;
                    catch
                        fprintf(fileID,'t=%1.2f p1=%1.2f p2=%1.2f p3=%1.2f p6=%1.2f p7=%1.2f\n',[t p1 p2 p3 p6 p7]);
                        continue;
                    end
                end
            end
        end
    end
end
fclose(fileID);


