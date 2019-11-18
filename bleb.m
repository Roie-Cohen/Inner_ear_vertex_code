function g = bleb(g)
% gets a semi-organized lattice and finishes the organiztion
global timer


g.fa = [1 1 1 1.2 1];
g = redistributeAreas(g);
g.bc = 0;

sid = 2;
vid = VideoWriter(strcat('../simulations/bleb_sim(', num2str(sid),')','.avi'));
open(vid);

timer = 0;
% parameters: area, tension, perimeter, mechanical feedback, HCs vdv, PCs
% pull, axial flow of SCs
g.paras = [3 ; 0.0; 0.15; 0; 0; 0; 0];
g = relaxLattice(g,0,vid);
T1prob = 1;
T1eps = 0.05;
asf = 0;
for t =1:30
    g = findTransitions(g,T1eps,T1prob,0.02); %0.05, 0.02, 0.02
    g = updateParameters(g);
    g = relaxLattice(g,5,vid); %50
    
    g = killSmallCells(g,0.1);
    g = updateParameters(g);
%     g = findUnstableBonds(g);
    g = relaxLattice(g,5,vid); %200
    
    disp(num2str(t))
end
close(vid);
end