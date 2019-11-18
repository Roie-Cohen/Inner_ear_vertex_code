function dE = dHflow2(g, c, R)
% simulates the flow of the outer region along the base-apex axis
    
dE = zeros(2*length(g.verts),1);

% if g.populations(c) == 2
    
   vidx = g.bonds(g.cells{c+1},1);
   vn = length(vidx); % number of vertices in cell c
   
   PCs = find(g.populations == 4);
   PCs_pos = g.centroid(PCs, :);
   PCs_y = mean(PCs_pos(:,2));
   
   cpos = g.centroid(c, :);
   
   D = cpos(2) - PCs_y;
   
   if D > 0 && D < 1.25*R
       if D < R
           for k=1:vn
               dE(2*vidx(k)-1) = dE(2*vidx(k)-1) - D/R;
           end
       end
       for k=1:vn
           dE(2*vidx(k)) = dE(2*vidx(k)) - 0*D/R;
       end
   end
   
% end
dE = g.paras(7)*dE;
end