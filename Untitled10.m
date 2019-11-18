global timer
wb = waitbar(0,'winter is coming');
for i=1:30
load(['stage1_complete-211118\lat(',num2str(i) ,')_step(30).mat'],'g');
g = stage2differentiate(g);

g.area_feedback = 0;
g.fa = [1 1 1 1.2 1];
g = redistributeAreas(g);

timer = 0;
% parameters: area, tension, perimeter, mechanical feedback, HCs vdv, PCs
% pull, axial flow of SCs
g.paras = [4 ; 0.4; 0.1; 0; 1; 0; 0];
T1prob = 1;
T1eps = 0.1;
for t =1:101
    g = findTransitions(g,T1eps,T1prob,0.02); %0.05, 0.02, 0.02
    g = updateParameters(g);
    g = relaxLattice(g,5); %50

    g = killSmallCells(g,0.30);
    g = updateParameters(g);
    g = relaxLattice(g,5); %200
    
    disp(num2str(t))
    
    if mod(t-1,10) == 0
        save(['stage2_complete-211118/lat(', num2str(i),')_step(', num2str((t-1)/10) ,').mat'],'g');
    end
end

waitbar(i/30);
end
close(wb)