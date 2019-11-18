function g= disorderedLattice(nrow,ncol,with_trans)

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
g.areas = zeros(length(g.cells)-1,1);
LatticePresentation(g,0);
if(g.bc==1)
g = periodicBC(g,nrow,ncol);
end
g.area_feedback = 0;
g = rescale(g);
%g.dead =zeros(length(g.cells)-1,1);
hex_area = cellarea(g,1);

% assigning arbitrary preffered areas and polarities 
area_fac = 0.4;
% for i=1:length(g.cells)-1
%     if(g.dead(i)==0)
%         g.areas(i) = (1+area_fac*2*(rand-0.5))*hex_area;     
%     end
% end
g.areas = (1+area_fac*3*(rand(length(g.areas),1)-0.5))*hex_area;

% random tensions
g.bonds(:,5) = 1 + rand(length(g.bonds),1);
g.populations = zeros(length(g.cells)-1,1);
g.cmpr = ones(length(g.cells)-1,1);
g.ctrc = ones(length(g.cells)-1,1);
g.transitionedBonds = [0, 0];

sid = 100;
vid = VideoWriter(strcat('../simulations/', 'dis_lat', '(', num2str(sid),')','.avi'));
open(vid);

timer = 0;
% parameters: area, tension, perimeter, mechanical feedback, HCs vdv
g.paras = [1 ; 0.15; 0; 0; 0];
g = relaxLattice(g,10,vid);
T1prob = 0.3;
T1eps = 0.05;
if(with_trans)
    for t =1:80,
        g = findTransitions(g,T1eps,T1prob,0.02); %0.05, 0.02, 0.02
        g = relaxLattice(g,10, vid); %50
        g = killSmallCells(g,0.05); %0.1
        g = relaxLattice(g,10, vid); %200
        disp(num2str(t))
        if (mod(t,8) == 0) && t<75
            g.bonds(:,5) = 1 + 2*rand(length(g.bonds),1);
            g.areas = (1+area_fac*2*(rand(length(g.areas),1)-0.5))*hex_area;
        end
        if t==60
            area_fac = 0.2;
            g.paras(3)=0.01;
        end
        if t==74
            g.bonds(:,5) = 1;
        end
    end
end
close(vid);
disp('finished');
figure(4),LatticePresentation(g,0);
