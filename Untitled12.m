

% give the function lattice (g) 
% do voronoin
% for each HC get the distance of the bonds from the centroid
% the minimal distance is the radius of the HC
% draw the circles on the lattice and think

HCi = find(g.populations==3 & ~g.dead);
x = g.centroid(HCi,1);
y = g.centroid(HCi,2);
x = [x; x+2*pi; x-2*pi];
y = [y;y;y];

% LatticePresentation(g,1);
% hold on
% voronoi(x,y)

[ve, cc] = voronoin([x, y]);

nc = length(HCi);
Rs = zeros(nc,1); % radii of the HCs
d0 = norm([x(1)-x(2),y(1)-y(2)]); % some lengthscale
for i = 1:nc
    hc = HCi(i); % index of the HC
    nvi = length(cc{i}); % number of vertices
    dmin = d0;
    for v=1:nvi
        pos1 = ve(cc{i}(v),:); % position of the vertices of the bond
        pos2 = ve(cc{i}(mod(v,nvi)+1),:);
        m = (pos2(2)-pos1(2))/(pos2(1)-pos1(1)); % bond slope
        n = pos1(2) - m*pos1(1); % crossing with x=0
        d = abs((m*x(i)-y(i)+n)/sqrt(m^2+1)); % distance of the cell centroid from the bond
        if d < dmin, dmin = d; end
    end
    Rs(i) = dmin;
end

LatticePresentation(g,0);
hold on
for i = 1:nc
    draw_circle(x(i),y(i),Rs(i));
end

%%
%HCs real area
HCs_real_area = pi*Rs.^2;

SCi = find(g.populations~=3 & ~g.dead);
SCs_real_area = zeros(length(SCi),1);
for i=1:length(SCi)
    sc = SCi(i);
    vidx = g.bonds(g.cells{sc+1},1);
    vert= getRelativePosition(g,vidx,sc); % vertices position
    xv = vert(:,1); yv = vert(:,2); % x and y coordinates
    xv = [xv ; xv(1)]; yv = [yv ; yv(1)];
    xp = min(xv) + (max(xv)-min(xv))*rand(1000,1); 
    yp = min(yv) + (max(yv)-min(yv))*rand(1000,1); 
    in = inpolygon(xp,yp,xv,yv);
    poly_area = sum(in);
    
    bonds = g.cells{sc+1};
    in_real = in; % logical array, checks if the point is in the SC but out of the circles
    for b=1:length(bonds)
        neigh = g.bonds(bonds(b),4);
        if g.populations(neigh)~= 3, continue; end
        ind = find(HCi == neigh,1);
        in_real = in_real & ( (xp-x(ind)).^2 + (yp-y(ind)).^2 > Rs(ind)^2 );
    end
    real_area = sum(in_real);
    
    SCs_real_area(i) = (real_area/poly_area)*cellarea(g, sc);
end

%%

LatticePresentation(g,0);
hold on
areas = zeros(length(g.cells)-1,1);
HCs = find(g.populations==3 & ~g.dead);
yhcs = g.centroid(HCs,2);
PC_h = mean(g.centroid(g.populations==4 & ~g.dead, 2));
for ii = -1:2:1 %inner and outer
    HCi = HCs(sign(yhcs-PC_h)==ii);
    x = g.centroid(HCi,1);
    y = g.centroid(HCi,2);
    x = [x; x+2*pi; x-2*pi];
    y = [y;y;y];
    
    [ve, cc] = voronoin([x, y]);
    
    nc = length(HCi);
    Rs = zeros(nc,1); % radii of the HCs
    d0 = norm([x(1)-x(2),y(1)-y(2)]); % some lengthscale
    for i = 1:nc
        nvi = length(cc{i}); % number of vertices
        dmin = d0;
        for v=1:nvi
            pos1 = ve(cc{i}(v),:); % position of the vertices of the bond
            pos2 = ve(cc{i}(mod(v,nvi)+1),:);
            m = (pos2(2)-pos1(2))/(pos2(1)-pos1(1)); % bond slope
            n = pos1(2) - m*pos1(1); % crossing with x=0
            d = abs((m*x(i)-y(i)+n)/sqrt(m^2+1)); % distance of the cell centroid from the bond
            if d < dmin, dmin = d; end
        end
        Rs(i) = dmin;
        draw_circle(x(i),y(i),Rs(i));
    end   
end