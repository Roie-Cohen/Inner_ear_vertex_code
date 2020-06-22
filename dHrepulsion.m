function dE = dHrepulsion(g,i)
    dE = zeros(2*length(g.verts),1);
    if g.populations(i) == 3 % if this is a HC
    % parameters for Lennard-Jones potential (k>l)
    kappa = 12; l = 0; r_min = 0.8; sigma = r_min*(l/kappa)^(1/(kappa-l)); 
    if l==0, sigma=0.7; end % if the potential is only repulsion.
    dEvdv = @(r) max(l*sigma^l./(r.^(l+1)) - kappa*sigma^kappa./(r.^(kappa+1)), -1000); % limiting the energy for the digit accuracy
    
        vidx = g.bonds(g.cells{i+1},1);
        vert = getRelativePosition(g,vidx);
        vn = length(vidx); % number of vertices in cell n
        Cn = g.centroid(i,:);
        % inner HCs don't experience the repulsion potential
        PCs = g.populations == 4;
        PCs_pos = g.centroid(PCs, :);
        PCs_y = mean(PCs_pos(:,2));
        if Cn(2) < PCs_y, return; end
        iHCs = find(g.populations == 3 & ~g.dead); % HCs indices
        iHCs(iHCs == i) = [];  % removing the current cell
        nHCs = length(iHCs);
        A = cellarea(g, i);
        for j = 1:nHCs
            m = iHCs(j);    % HC index
            Cm = g.centroid(m, :);
            if Cm(2) < PCs_y, continue; end % inner HCs don't play
            r = Cn - Cm;
%             if norm(r)>sigma, continue; end % only steric repultion
            % if we're with a periodic boundary conditions we need to check
            % if the distance between 'n' and 'm' is actually closer
            if g.bc == 1
                % the scaling of the lattice in periodic BC is 2pi, meaning
                % the height and width of the lattice is 2pi. 
                % if there's a closer matching vertex in the periodic
                % lattice, then r is replaced with the vector pointing from
                % the closer point.
                r(1) = r(1)-(2*pi-abs(r(1)) < abs(r(1)))*sign(r(1))*2*pi;
                r(2) = r(2)-(2*pi-abs(r(2)) < abs(r(2)))*sign(r(2))*2*pi;
            end
            if norm(r) > pi
                continue;
            end
            % dEnm = dE/dr * dr/dx
            dr = r/norm(r);
            dEnmdr = dEvdv(norm(r));
            for k=1:vn
                kn = mod(k,vn)+1; % next vertex
                kp = mod(k-2,vn)+1; % previous vertex
                r0 = vert(kp,:); r1 = vert(k,:); r2 = vert(kn,:); % positions of vertices k-1, k, k+1 respectively
                xp = r0(1); x = r1(1); xn = r2(1); yp = r0(2); y = r1(2); yn = r2(2);
                dAdx = 0.5*( yn-yp );
                dAdy = 0.5*( xp-xn );
                dXcmdx = ( -dAdx*Cn(1) + (1/6)*(2*x*(yn-yp)+y*(xp-xn)+xn*yn-xp*yp) )/A;
                dYcmdx = ( -dAdx*Cn(2) + (1/6)*(yn*(y+yn)-yp*(y+yp)) )/A;
                dXcmdy = ( -dAdy*Cn(1) + (1/6)*(xp*(x+xp)-xn*(x+xn)) )/A;
                dYcmdy = ( -dAdy*Cn(2) + (1/6)*(2*y*(xp-xn)+x*(yn-yp)+yp*xp-yn*xn) )/A;
                dEnm = dEnmdr*[dr(1)*dXcmdx+dr(2)*dYcmdx, dr(1)*dXcmdy+dr(2)*dYcmdy];
                dE(2*vidx(k)-1) = dE(2*vidx(k)-1) - dEnm(1);
                dE(2*vidx(k)) = dE(2*vidx(k)) - dEnm(2);
            end
        end
        dE = g.paras(4)*dE;

    end
end
