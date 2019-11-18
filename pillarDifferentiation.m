function g = pillarDifferentiation(g, R)
% gets the distance of differentiation from the pillar cells R
% cells in a distance of R from the pillar cells will differentiate into
% HCs.

if isfield(g,'populations')
   PCs = find(g.populations == 4);
   PCs_pos = cellCenter(g, PCs);
   nc = length(g.cells)-1;
   ci = randperm(nc);
   ci = ci(g.populations(ci) == 2 & ~g.dead(ci)); 
   for i=1:length(ci)
       % check if there's an adjecent HC
       if ~isempty(find(g.populations(g.bonds(g.cells{ci(i)+1},4))==3,1))
           continue;
       end
       % check the distance from the PCs
       cpos = cellCenter(g, ci(i));
       Ds = sqrt(sum((PCs_pos - cpos).^2,2));
       [m, mi] = min(Ds);
       dy = cpos(2) - PCs_pos(mi, 2);
       if dy < R && dy > 0
           g.populations(ci(i)) = 3;     
       end
   end
    
end

end