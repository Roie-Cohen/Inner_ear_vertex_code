function dE=dHperimeter(g,i)
dE = zeros(2*length(g.verts),1);
vidx=g.bonds(g.cells{i+1},1); % an array of the vertices indices of the cell
vert = getRelativePosition(g,vidx,i); % the position of the vertices
nb=length(vidx);
epsilon_g = 0; %global polarity parameter
epsilon_l = 0; %local polarity parameter
for j=1:nb
  prev = mod(j-2,nb)+1; % the previous vertex index (in vert)
  next = mod(j,nb)+1;   % the next vertex
  prev_theta = atan( (vert(j,2)-vert(prev,2))/(vert(j,1)-vert(prev,1)) ); 
  next_theta = atan( (vert(j,2)-vert(next,2))/(vert(j,1)-vert(next,1)) );
 n1 = norm(vert(j,:)-vert(prev,:));
 n2 = norm(vert(j,:)-vert(next,:));
 if(n1>0.0001),
     border_cells = g.bonds(g.cells{i+1}(prev), [3, 4]); % the cells separated by the bond defined by the previous vertex
     local_polarity = 1 + epsilon_l*( cos( g.polarity(border_cells(1)) - prev_theta )^2 );
     if border_cells(2) ~= 0,
        local_polarity = local_polarity + epsilon_l*( cos( g.polarity(border_cells(2)) - prev_theta )^2 );
     end
     global_polarity = 1+epsilon_g*cos(prev_theta)^2;
 dE(2*vidx(j)-1) = local_polarity*global_polarity*(vert(j,1)-vert(prev,1))/n1;
 dE(2*vidx(j)) =   local_polarity*global_polarity*(vert(j,2)-vert(prev,2))/n1;
 else
   dE(2*vidx(j)-1) = 1;
 dE(2*vidx(j)) = 1;
 end
 
  if(n2>0.0001),
      border_cells = g.bonds(g.cells{i+1}(j), [3, 4]); % the cells separated by the bond defined by the next vertex
      local_polarity = 1 + epsilon_l*( cos( g.polarity(border_cells(1)) - next_theta )^2 );
      if border_cells(2) ~= 0,
          local_polarity = local_polarity + epsilon_l*( cos( g.polarity(border_cells(2)) - next_theta )^2 );
      end
      global_polarity = 1+epsilon_g*cos(next_theta)^2;
      dE(2*vidx(j)-1) =  dE(2*vidx(j)-1)+local_polarity*global_polarity*(vert(j,1)-vert(next,1))/n2;
      dE(2*vidx(j)) =   dE(2*vidx(j))+ local_polarity*global_polarity*(vert(j,2)-vert(next,2))/n2;
 else
   dE(2*vidx(j)-1) =  dE(2*vidx(j)-1)+1;
 dE(2*vidx(j)) = dE(2*vidx(j))+1;
 end
 end
dE = dE* g.paras(2);
end