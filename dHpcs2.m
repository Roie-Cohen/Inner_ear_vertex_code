function dE = dHpcs2(g, c)
% gets the attraction energy differential for cell c.
% E(c) = d^2 where d is the distance from cell c to the line on PCs.
dE = zeros(2*length(g.verts),1);

if g.populations(c) == 3
%     if ~isempty(find(g.populations(g.bonds(g.cells{c+1},4))== 4,1)), return; end
   vidx = g.bonds(g.cells{c+1},1);
   vn = length(vidx); % number of vertices in cell c
   
   PCs = find(g.populations == 4);
   PCs_pos = g.centroid(PCs, :);
   Lin_fit = polyfit(PCs_pos(:,1), PCs_pos(:,2), 1); % line defining the pillar cells axis
   cpos = g.centroid(c,:);
   m = Lin_fit(1); n = Lin_fit(2); x = cpos(1); y = cpos(2);
   
   D = (y-m*x-n)/sqrt(m^2+1);
   
   vmax = [0 0];
   kmax = 0;
   for k=1:vn
       % finding the angle of the bond going out from vetex vids(k) (away
       % from the cell). The angle relative to +x axis.
       v1 = vidx(k);
       nbs = find(g.bonds(:,1)==v1);
       nvs = g.bonds(nbs, 2);
       v2 = nvs(~ismember(nvs, vidx));
       dv = g.verts(v2,1:2)-g.verts(v1,1:2);
       if norm(dv) < 0.01, continue; end
       vhat = dv/norm(dv);
       if vhat(2)>=0, continue; end
       ncs = find_bond_edge_cells(g,nbs(nvs==v2));
       ncp = g.populations(ncs(ncs~=c));
       if ncp == 3 && norm(dv)<0.1, continue; end
       dE(2*vidx(k)-1) = dE(2*vidx(k)-1) + (2*D+1)*vhat(1)*vhat(2);
       dE(2*vidx(k)) = dE(2*vidx(k)) + (2*D+1)*vhat(2)^2;
   end


end
dE = g.paras(6)*dE;
end