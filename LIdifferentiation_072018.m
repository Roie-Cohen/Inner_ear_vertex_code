function g = LIdifferentiation(g, R)
% gets the distance of differentiation from the pillar cells R
% cells in a distance of R from the pillar cells will differentiate into
% HCs.

single_diff_prop = 1; % probability for differentiation of individual cells (0.75 works good)
if isfield(g,'populations')
   PCs = find(g.populations == 4);
   PCs_pos = cellCenter(g, PCs);
   nc = length(g.cells)-1;
   ci = randperm(nc);
   ci = ci(g.populations(ci) == 2 & ~g.dead(ci)); 
   snHCs = zeros(length(ci), 1);    % second neighboring HCs
   for i=1:length(ci)
       % check if there's an adjecent HC
       nbrs = g.bonds(g.cells{ci(i)+1},4);
       if ~isempty(find(g.populations(nbrs)==3,1))
           continue;
       end
       % check the distance from the PCs
       cpos = cellCenter(g, ci(i));
       Ds = sqrt(sum((PCs_pos - cpos).^2,2));
       [m, mi] = min(Ds);
       dy = cpos(2) - PCs_pos(mi, 2);
       if dy < R && dy > 0
           % finds how much HCs are second neighbors to the cell
           snbrs = [];
           for j=1:length(nbrs)
               snbrs = [snbrs; g.bonds(g.cells{nbrs(j)+1}, 4)];
           end
           snbrs = snbrs(snbrs~=ci(i));
           snbrs = unique(snbrs);
           snHCs(ci(i)) = sum(g.populations(snbrs)==3);
           
       end
   end
   mx = max(snHCs);
   if mx~=0 && rand() < single_diff_prop
       mxs = find(snHCs == mx);
       g.populations(mxs(randi(length(mxs)))) = 3;
   end
    
end

end