function g = LIdifferentiation(g, R)
% gets the distance of differentiation from the pillar cells R
% cells in a distance of R from the pillar cells will differentiate into
% HCs.

if sum(g.populations == 3) >= 27, return; end

flag = 0;
single_diff_prop = 1; % probability for differentiation of individual cells (0.75 works good)
if isfield(g,'populations')
%    PCs = find(g.populations == 4);
%    PCs_pos = cellCOM(g, PCs); 
   PCs_pos = g.centroid(g.populations == 4, :);
   nc = length(g.cells)-1;
   ci = randperm(nc);
   ci = ci(g.populations(ci) == 2 & ~g.dead(ci)); 
%    snHCs = zeros(length(ci), 1);    % second neighboring HCs
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
           flag = 1;
       end
   end
    
end

if flag, g = redistributeAreas(g); end

end