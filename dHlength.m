
function dE=dHlength(g,i)
% calculates the energy deferential derived from the bonds length and
% perimeter elements in the energy function.
l_fac = g.paras(2);
L_fac = g.paras(3)*g.ctrc(i);
dE = zeros(2*length(g.verts),1);
L = cellPerimeter(g,i);
vidx=g.bonds(g.cells{i+1},1); % an array of the vertices indices of the cell
vert = getRelativePosition(g,vidx,i); % the position of the vertices
nb=length(vidx);
epsilon_g = 0; %global polarity parameter
for j=1:nb
    prev = mod(j-2,nb)+1; % the previous vertex index (in vert)
    next = mod(j,nb)+1;   % the next vertex
    prev_theta = atan( (vert(j,2)-vert(prev,2))/(vert(j,1)-vert(prev,1)) );
    next_theta = atan( (vert(j,2)-vert(next,2))/(vert(j,1)-vert(next,1)) );
    n1 = norm(vert(j,:)-vert(prev,:));
    n2 = norm(vert(j,:)-vert(next,:));
    if(n1>0.0001),
        border_cells = g.bonds(g.cells{i+1}(prev), [3, 4]); % the cells separated by the bond defined by the previous vertex
        global_polarity = 1+epsilon_g*sin(prev_theta)^2;
        popfac = g.bonds(g.cells{i+1}(prev), 5);
        dE(2*vidx(j)-1) = (L_fac*L + l_fac*popfac*global_polarity)*(vert(j,1)-vert(prev,1))/n1;
        dE(2*vidx(j)) =   (L_fac*L + l_fac*popfac*global_polarity)*(vert(j,2)-vert(prev,2))/n1;
    else
        dE(2*vidx(j)-1) = 1;
        dE(2*vidx(j)) = 1;
    end
    
    if(n2>0.0001),
        border_cells = g.bonds(g.cells{i+1}(j), [3, 4]); % the cells separated by the bond defined by the next vertex
        global_polarity = 1+epsilon_g*sin(next_theta)^2;
        popfac = g.bonds(g.cells{i+1}(j), 5);
        dE(2*vidx(j)-1) =  dE(2*vidx(j)-1) + (L_fac*L + l_fac*popfac*global_polarity)*(vert(j,1)-vert(next,1))/n2;
        dE(2*vidx(j)) =   dE(2*vidx(j)) + (L_fac*L + l_fac*popfac*global_polarity)*(vert(j,2)-vert(next,2))/n2;
    else
        dE(2*vidx(j)-1) =  dE(2*vidx(j)-1)+1;
        dE(2*vidx(j)) = dE(2*vidx(j))+1;
    end
end
% dE = dE* g.paras(2);
end
