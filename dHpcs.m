function dE = dHpcs(g, c)
% gets the attraction energy differential for cell c.
% E(c) = d^2 where d is the distance from cell c to the line on PCs.
dE = zeros(2*length(g.verts),1);

if g.populations(c) == 3
%     if ~isempty(find(g.populations(g.bonds(g.cells{c+1},4))== 4,1)), return; end
   vidx = g.bonds(g.cells{c+1},1);
   vn = length(vidx); % number of vertices in cell c
   
   PCs = g.populations == 4;
   PCs_pos = g.centroid(PCs, :);
   Lin_fit = polyfit(PCs_pos(:,1), PCs_pos(:,2), 1); % line defining the pillar cells axis
   cpos = g.centroid(c,:);
   m = Lin_fit(1); n = Lin_fit(2); x = cpos(1); y = cpos(2);
   
   D = (y-m*x-n)/sqrt(m^2+1);
   Dx = -2*(D)*m/(sqrt(m^2+1)*vn);
   Dy = 2*(D)/(sqrt(m^2+1)*vn);
   
   for k=1:vn
       dE(2*vidx(k)-1) = dE(2*vidx(k)-1) + Dx;
       dE(2*vidx(k)) = dE(2*vidx(k)) + Dy;
   end
 
%    lead_v = vidx(g.verts(vidx,2) <= cpos(2));
%    ln = length(lead_v);
%    Dy = Dy*vn/ln;
%    for k=1:ln
%        dE(2*lead_v(k)) = dE(2*lead_v(k)) + Dy;
%    end

end
dE = g.paras(6)*dE;
end