function g= ImportedLattice(with_trans)

g = getGfromCIA();
g = rescale(g);
g.area_feedback = 1;

% assigning arbitrary preffered areas and polarities 
area_fac = 0;
for i=1:length(g.cells)-1,
    if(g.dead(i)==0)
%         g.areas(i) = (1+area_fac*2*(rand-0.5))*cellarea(g,i);
        g.polarity(i) = (2*rand-1)*pi;
    end
end

%area-feedback model
if (g.area_feedback==1)
    g.signal = zeros(length(g.cells)-1,1);
    g.delta = zeros(length(g.cells)-1,1);
    g = initiateDeltaSignal(g);
end

figure(1), LatticePresentation(g,1);

switch with_trans
    case 0, pre = 'r';
    case 1, pre = 'i';
end
switch g.bc
    case 0, bc='np';
    case 1, bc='p';
    case 2, bc='st';
end
sid = 10;
vid = VideoWriter(strcat('../simulations/', pre, '_imp', bc,'(', num2str(sid),')','.avi'));
open(vid);

g.paras = [1 ; 0.7; 0.1; 1.5];
g = relaxLattice(g,500,vid);
if(with_trans)
    iter = 6;
    t1prob = 0.01;
    for t =1:iter
        g = findTransitions(g,0.2,t1prob,0.02); %1, 0.02, 0.02
        g = relaxLattice(g,50,vid); %50
        g = killSmallCells(g,0.1,2);
        g = findUnstableBonds(g);
        g = relaxLattice(g,100,vid); %200   
        disp(t)
    end
    g = relaxLattice(g,200,vid);
end
close(vid);
figure(3),LatticePresentation(g,0);
