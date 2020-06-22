function g = LIdifferentiation(g, R)
% gets the distance of differentiation from the pillar cells R
% cells in a distance of R from the pillar cells will differentiate into
% HCs.

% limiting the simulation to 21 OHCs + 6 IHCs
if sum(g.populations == 3) >= 27, return; end 

PCs_pos = g.centroid(g.populations == 4, :);
nc = length(g.cells)-1;
ci = randperm(nc);
ci = ci(g.populations(ci) == 2 & ~g.dead(ci));
for i=1:length(ci)
    % check if there's an adjecent HC
    nbrs = g.bonds(g.cells{ci(i)+1},4);
    if ~isempty(find(g.populations(nbrs)==3,1))
        if surroundingHCsFrac(g, ci(i)) > 0.05
            continue;
        end
    end
    % check the distance from the PCs
    cpos = g.centroid(ci(i), :); %cellCOM(g, ci(i));
    Ds = sqrt(sum((PCs_pos - cpos).^2,2));
    [~, mi] = min(Ds);
    dy = cpos(2) - PCs_pos(mi, 2);
    if dy < R && dy > 0
        g.populations(ci(i)) = 3;
    end
end

end