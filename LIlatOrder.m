function g = LIlatOrder(sid)
% generates an irrgelular lattice  with nrow rows and ncol columns.
% If with_trans is false, all cells have exactly  6 neighbors,
% otherwise some cells may be removed.
warning off
global ctypes timer T1delay noise
ctypes = [];
T1delay = 0;

% load(['lattices/LIlat_autodiff', num2str(sid)],'g');
load(['lattices/LIlat_autodiff1'],'g');
% load('lattices\main_lat271118-after(7).mat','g');

g.transitionedBonds = [0, 0];
g.area_feedback = 0;
g.fa = [1 1 0.7 1.2 1 1];  % Hensen:SC:HCs:pillar area ratio
g = redistributeAreas(g);
g.cmpr = ones(length(g.cells) - 1, 1); % compressibility of cells
g.ctrc = ones(length(g.cells) - 1, 1); % contractility of cells
g.area_pert = zeros(length(g.cells) - 1, 1); % logial array, if a cell's area is perturbed
g.centroid = zeros(length(g.cells) - 1, 2); % centroids of cells

%area-feedback model
if (g.area_feedback==1)
    g.signal = zeros(length(g.cells)-1,1);
    g.delta = zeros(length(g.cells)-1,1);
    g = initiateDeltaSignal(g);
end

vid = VideoWriter(strcat('../simulations/stage1-',num2str(sid),'.avi'));
vid.Quality = 100;
open(vid);

timer = 0;
% parameters: area, tension, perimeter, mechanical feedback, HCs vdv, PCs
% pull, axial flow
g.paras = [6 ; 0.3; 0.15; 0; 1; 0.5; 0.4]; %[9 ; 0.3; 0.15; 0; 1; 0.5; 0.4]; 
g = relaxLattice(g,0,vid);
T1prob = 1; % 0.5
T1eps = 0.1;
g = LIdifferentiation(g, 2.4);
noise = rand(1, 2*length(g.verts));
for t =1:200
    
    g = findTransitions(g,T1eps,T1prob,0.02); %0.05, 0.02, 0.02
    g = updateParameters(g);
    g = relaxLattice(g,5,vid); %50
    
    g = killSmallCells(g,0.1);
    g = updateParameters(g);
%     g = findUnstableBonds(g);
    g = relaxLattice(g,5,vid); %200
    
    g = LIdifferentiation(g, 2.6);
    g = updateParameters(g);
    g = relaxLattice(g,10,vid); %50
    
    disp(num2str(t))
end


close(vid);
figure(4),LatticePresentation(g,0);
