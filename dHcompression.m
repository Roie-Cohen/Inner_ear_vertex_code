function dE = dHcompression(g, c)
% gets the attraction energy differential for cell c.
% E(c) = d^2 where d is the distance from cell c to the line on PCs.
dE = zeros(2*length(g.verts),1);

if g.populations(c) == 3
    %     if ~isempty(find(g.populations(g.bonds(g.cells{c+1},4))== 4,1)), return; end
    vidx = g.bonds(g.cells{c+1},1);
    vn = length(vidx); % number of vertices in cell c
    vert = getRelativePosition(g,vidx);
    PCs = g.populations == 4;
    PCs_pos = g.centroid(PCs, :);
    PCs_y = mean(PCs_pos(:,2));
    cpos = g.centroid(c,:);
    if cpos(2) < PCs_y, return; end % force don't act on inner HCs
    A = cellarea(g, c);
    xc = cpos(1); yc = cpos(2);
    
    D = abs(yc-PCs_y);
    
    for k=1:vn
        kn = mod(k,vn)+1; % next vertex
        kp = mod(k-2,vn)+1; % previous vertex
        r0 = vert(kp,:); r1 = vert(k,:); r2 = vert(kn,:); % positions of vertices k-1, k, k+1 respectively
        xp = r0(1); x = r1(1); xn = r2(1); yp = r0(2); y = r1(2); yn = r2(2);
        dAdx = 0.5*( yn-yp );
        dAdy = 0.5*( xp-xn );
        dYcmdx = ( -dAdx*yc + (1/6)*(yn*(y+yn)-yp*(y+yp)) )/A;
        dYcmdy = ( -dAdy*yc + (1/6)*(2*y*(xp-xn)+x*(yn-yp)+yp*xp-yn*xn) )/A;
        dE(2*vidx(k)-1) = dE(2*vidx(k)-1) - 2*D*dYcmdx;
        dE(2*vidx(k)) = dE(2*vidx(k)) - 2*D*dYcmdy;
    end
    
end
dE = g.paras(5)*dE;
end