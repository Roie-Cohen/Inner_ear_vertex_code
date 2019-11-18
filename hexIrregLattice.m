function g= hexIrregLattice(nrow,ncol,with_trans)
% generates an irrgelular lattice  with nrow rows and ncol columns.
% If with_trans is false, all cells have exactly  6 neighbors,
% otherwise some cells may be removed.
global ctypes timer T1delay
ctypes = [];
T1delay = 0;


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



disp('Number of cells: ');disp(length(cc));
g(1) = GLattConversion(cc,vertices);
g.xboundary=zeros(length(g.cells)-1,2);
g.yboundary=zeros(length(g.cells)-1,2);
g.bc=1;
g.scale =eye(2);
g.dead =zeros(length(g.cells)-1,1);
LatticePresentation(g,0);
if(g.bc==1)
g = periodicBC(g,nrow,ncol);
end

g = rescale(g);
g.dead =zeros(length(g.cells)-1,1);
hex_area = cellarea(g,1);

disordered = 1;
if (disordered)
    load('lattices/dis_lat2');
%     load('lattices/LIlat1');
end

g.area_feedback = 0;
g.fa = [1 1 1 1.2];  % Hensen:SC:HCs:pillar area ratio

% assigning arbitrary preffered areas and polarities 
area_fac = 0;
for i=1:length(g.cells)-1,
    if(g.dead(i)==0)
        %g.areas(i) = (1+area_fac*2*(rand-0.5))*cellarea(g,i);
        g.areas(i) = (1+area_fac*2*(rand-0.5))*hex_area;
%         g.areas(i) = (1-0.6*( rand()<0.8 ))*hex_area;
        g.polarity(i) = (2*rand-1)*pi;
        
    end
end

populations = 1;
if populations == 1
    g = differentiate(g, 4);
    g = redistributeAreas(g);
end

        
% stretching the tissue
strch_amount = 1;
if g.bc == 2,
    g = stretchTissue(g, strch_amount);
end

%area-feedback model
if (g.area_feedback==1)
    g.signal = zeros(length(g.cells)-1,1);
    g.delta = zeros(length(g.cells)-1,1);
    g = initiateDeltaSignal(g);
end

figure(3), LatticePresentation(g,1);

switch with_trans
    case 0, pre = 'r';
    case 1, pre = 'i';
end
switch g.bc
    case 0, bc='np';
    case 1, bc='p';
    case 2, bc='st';
end
if (disordered),
    lattype = 'dis';
else
    lattype = 'h12x12';
end
sid = 133;
vid = VideoWriter(strcat('../simulations/', pre, lattype, bc,'(', num2str(sid),')','.avi'));
open(vid);

timer = 0;
% parameters: area, tension, perimeter, mechanical feedback, HCs vdv
g.paras = [6 ; 0.1; 0.1; 0; 5];
g = relaxLattice(g,0,vid);
T1prob = 0.50;
T1eps = 0.05;
asf = 0;
if(with_trans)
    for t =1:200,
        g = findTransitions(g,T1eps,T1prob,0.02); %0.05, 0.02, 0.02
        g = relaxLattice(g,10,vid); %50
        g = killSmallCells(g,0.1);
        g = findUnstableBonds(g);
        g = relaxLattice(g,10,vid); %200
        disp(num2str(t))
%         if t<=50 
%             g.fa(3) = 1 + (0.5/50)*t;
%             g = redistributeAreas(g);
%         end
        
    end

%     g = relaxLattice(g,200,vid);
end
close(vid);
figure(4),LatticePresentation(g,0);
