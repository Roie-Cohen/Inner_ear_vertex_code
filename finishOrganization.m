function g = finishOrganization(g)
% gets a semi-organized lattice and finishes the organiztion
global timer

g = stage2differentiate(g);

%area-feedback model
g.area_feedback = 1;
if (g.area_feedback==1)
    g.signal = zeros(length(g.cells)-1,1);
    g.delta = zeros(length(g.cells)-1,1);
    g = initiateDeltaSignal(g);
end
g.fa = [1 1 1 1.2 1];
g = redistributeAreas(g);

sid = 8;
vid = VideoWriter(strcat('../simulations/semi_org_lat(', num2str(sid),')','.avi'));
vid.Quality = 100;
open(vid);

timer = 0;
% parameters: area, tension, perimeter, mechanical feedback, HCs vdv, PCs
% pull, axial flow of SCs
g.paras = [4 ; 0.4; 0.1; 0; 1; 0; 0];
g = relaxLattice(g,0,vid);
T1prob = 1;
T1eps = 0.1;
asf = 0;
for t =1:100
    g = findTransitions(g,T1eps,T1prob,0.02); %0.05, 0.02, 0.02
    g = updateParameters(g);
    g = relaxLattice(g,5,vid); %50
    
%     if mod(t,10) == 0
%         answer = inputdlg({'cell:'});
%         ccc = str2num(answer{1});
%         if ccc ~= 0
%             g.areas(ccc) = 100;
%         end
%         save(['lattices/271118-t=',num2str(t) ,'.mat'], 'g');
%     end
    g = killSmallCells(g,0.20);
    g = updateParameters(g);
    g = relaxLattice(g,5,vid); %200
    
    disp(num2str(t))
end
close(vid);
end