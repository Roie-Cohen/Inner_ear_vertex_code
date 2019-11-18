function g = SignalAreaFeedback(g)
if(g.area_feedback)
    
    % defining the parameters
    p = struct('tau1',10,'tau2',0.05,'gamma',0.1,'h',4,'k1',3,'k2',10,'k3',15,'wm',0.01,'wr',10,'wh',50,'ws',30);
    
    % calculating the signal for each cell
    g.signal = zeros(length(g.cells)-1,1);
    cc = 1:length(g.cells)-1;
    for k=cc
        bonds = g.cells{k+1};
        vidx = g.bonds(bonds,1);                % an array of the vertices indices of the cell
        vert = getRelativePosition(g,vidx,k);   % the position of the vertices
        nb = length(bonds);    % number of bonds
        for j=1:nb
            neigh_cell = g.bonds(bonds(j),4);   % neighbor cell
            if neigh_cell~=0 && g.populations(neigh_cell)==3
                nextv = mod(j,nb)+1;                % the next vertex index
                lij = norm(vert(j,:)-vert(nextv,:));
                pj = cellPerimeter(g, neigh_cell);
                if pj<0.01, continue; end
                g.signal(k) = g.signal(k) + (lij/pj)*g.delta(neigh_cell);
            end
        end
        coeff = p.k2^p.h/(p.k2^p.h + g.delta(k)^p.h);
        g.signal(k) = g.signal(k)*coeff;
    end
    
    f = (p.k1^p.h)./(p.k1^p.h + g.signal.^p.h);
    
    eps = g.paras(4)*0.01;
    
    % advancing delta
    dD = p.tau1.*f - p.gamma.*g.delta;
    g.delta = g.delta + eps*dD;
    
    % advancing vertices positions
    nv = length(g.verts(:,1));
    dv = zeros(nv,2);
    for k=1:nv
        b = find(g.bonds(:,1)==k);
        vl = g.bonds(b ,2);
        for l=1:length(vl)
            sc = g.bonds(b(l),3:4);
%             if ~ismember(0,sc)
                if ismember(0,sc) || ismember(1,g.populations(sc)) 
                    dv(k,:) = [0 0];
                    break;
                end
                if sum(ismember([2 3],g.populations(sc))) == 2
                    c = max(g.signal(sc));
                    w = p.wm + p.wr*p.k3^p.h/(p.k3^p.h + c^p.h);
                else
                    if sum(g.populations(sc)==[3;3])==2
                        w = p.wh;
                    else
                        if (sum(g.populations(sc)==[2;2])==2) 
                            w = p.ws;
                        else
                            w = 0;
                        end
                    end
                end
                vpos = getRelativePosition(g,[k vl(l)]);
                dv(k,:) = dv(k,:) - w*( vpos(1,:) - vpos(2,:) );
%             end
        end
    end
    dv = p.tau2.*dv;
    g.verts(:,1:2) = g.verts(:,1:2) + eps.*dv;
    
end
end