% generate random lattices
warning off
global timer T1delay noise
clear g
nrow = 12;
ncol = 12;
[cc,ve] = mHexLattice(nrow+2,ncol+2,0);
ve(~isfinite(ve))=nan;
 
% remove unused vertices
vidx = zeros(length(ve),1);
for i = 1:length(cc)
    vidx(cc{i}) = 1;
end
vindex = find(vidx == 1);
vertices = ve(vindex,1:2);
vshift = zeros(length(ve),1);
vshift(vindex) = 1:length(vindex);
for i = 1:length(cc)
    cc{i} = vshift(cc{i})';
    cc{i}(end+1) = cc{i}(1);
end

g(1) = GLattConversion(cc,vertices);
g.xboundary=zeros(length(g.cells)-1,2);
g.yboundary=zeros(length(g.cells)-1,2);
g.bc=1;
g.scale =eye(2);
g.dead =zeros(length(g.cells)-1,1);
g.areas = zeros(length(g.cells)-1,1);
if(g.bc==1), g = periodicBC(g,nrow,ncol); end
g.area_feedback = 0;
g = rescale(g);
g.populations = zeros(length(g.cells)-1,1);
hex_area = cellarea(g,1);

g.transitionedBonds = [0, 0];
g.area_feedback = 0;
g.fa = [1 1 1 1 1 1];  % Hensen:SC:HCs:pillar area ratio
g = redistributeAreas(g);
g.cmpr = ones(length(g.cells) - 1, 1); % compressibility of cells
g.ctrc = ones(length(g.cells) - 1, 1); % contractility of cells
g.centroid = zeros(length(g.cells) - 1, 2); % centroids of cells
g.paras = [1 ; 0.1; 0; 0; 0; 0; 0];
g1 = g;


T1eps = 0.1;
wb = waitbar(0, 'Winter is coming...');

for l_ind=7:100
    
    try
        g = g1;
        T1delay = 0;
        timer = 0;
        noise = rand(1, 2*length(g.verts));  
        g.bonds(:,5) = 1 + rand(length(g.bonds),1); % randomize tensions
        g.areas = (1+2*(rand(length(g.areas),1)-0.5))*hex_area; % randomize areas
        T1prob = 1;
        for t =1:25
            g = findTransitions(g,T1eps,T1prob,0.02); %0.05, 0.02, 0.02
            g = relaxLattice(g,10);
            
            g = killSmallCells(g,0.1);
            g = relaxLattice(g,10);
            
            if mod(t,5) == 0
                g.bonds(:,5) = 1 + rand(length(g.bonds),1); % randomize tensions
                g.areas = (1+1*(rand(length(g.areas),1)-0.5))*hex_area; % randomize areas
            end
            if t==20
                g.areas(:) = hex_area;
                g.bonds(:,5) = 1;
                g.paras(3) = 0.05;
                g.bonds(:,5) = 1;
            end
        end
        fname = ['random_lattices-211118/rand_lat(', num2str(l_ind),').mat'];
        save(fname, 'g');
        waitbar(l_ind/100);
    catch
        disp(['error in random lattice #', num2str(l_ind)]);
        continue; 
    end
    
end
close(wb)

