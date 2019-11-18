function g = findUnstableBonds(g)
    % Searching for vertices that neighbors only HCs. 
    % the bond that touches the smallest SC is transitioning (T1)
    for i=1:length(g.verts(:,1))
        bidx = find(g.bonds(:,1) == i);
        neig_c = unique([g.bonds(bidx,3); g.bonds(bidx,4)]);
        if ismember(0, neig_c) || isempty(bidx)
            continue;
        end
        if ~ismember(0,ismember(g.populations(neig_c),3))
            neig_v = g.bonds(bidx,2);
            areas = zeros(length(neig_v),1);
            for j=1:length(neig_v)
                nbidx = find(g.bonds(:,2) == neig_v(j));
                nnc = unique([g.bonds(nbidx,3); g.bonds(nbidx,4)]);
                areas(j) = cellarea(g, nnc(~ismember(nnc,neig_c)));
            end
            [ord, id] = sort(areas);
            g = T1transition(g,bidx(id(1)));
        end
    end
    
end