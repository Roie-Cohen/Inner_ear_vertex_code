function dE=dHvdv(g,i)
if isfield(g,'populations')
    dE = zeros(2*length(g.verts),1);
    % parameters for Lennard-Jones potential (k>l)
    k = 12; l = 0; r_min = 0.8; sigma = r_min*(l/k)^(1/(k-l)); %r_min=0.8
    if l==0, sigma=0.7; end % if the potential is only repulsion.
    dEvdv = @(r) max(l*sigma^l./(r.^(l+1)) - k*sigma^k./(r.^(k+1)), -1000); % limiting the energy for the digit accuracy
    if g.populations(i) == 3 % if this is a HC
        vidx = g.bonds(g.cells{i+1},1);
        vn = length(vidx); % number of vertices in cell n
        Cn = g.centroid(i,:);
        iHCs = find(g.populations == 3 & ~g.dead); % HCs indices
        iHCs(iHCs == i) = [];  % removing the current cell
        nHCs = length(iHCs);
        for j = 1:nHCs
            m = iHCs(j);    % HC index
            Cm = g.centroid(m, :);
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
            dEnm = dEvdv(norm(r))*(r/norm(r))/vn;
            for k=1:vn
                dE(2*vidx(k)-1) = dE(2*vidx(k)-1) + dEnm(1);
                dE(2*vidx(k)) = dE(2*vidx(k)) + dEnm(2);
            end
        end
        dE = g.paras(5)*dE;

        % if the cell is touching a pillar cell then its vdv amplitude is
        % reduced
        if ~isempty(find(g.populations(g.bonds(g.cells{i+1},4))==4,1))
%             dE = 0.5*dE; 
        end
    end
end
end