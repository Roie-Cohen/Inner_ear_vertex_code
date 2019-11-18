function dE = dHflow3(g, c, R)
% simulates the flow of the HCs along the base-apex axis. 
% only the vertices in the front of the cell (in the direction of movement) is affected.

dE = zeros(2*length(g.verts),1);

if g.populations(c) == 3
    
   vidx = g.bonds(g.cells{c+1},1);
   vn = length(vidx); % number of vertices in cell c
   
   PCs = find(g.populations == 4);
   PCs_pos = cellCenter(g, PCs);
   PCs_y = mean(PCs_pos(:,2));
   
   cpos = cellCenter(g, c);
   
   D = cpos(2) - PCs_y;
   
   if D > 0 && D < 1.25*R
       if D < R
           lead_v = vidx(g.verts(vidx,1) >= cpos(1));
           ln = length(lead_v);
           for k=1:ln
               dE(2*lead_v(k)-1) = dE(2*lead_v(k)-1) - (D/R)/ln;
           end
       end
   end
   
 
end
dE = g.paras(7)*dE;
end