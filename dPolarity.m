
function dp=dPolarity(g,i)
dp = 0;
vidx=g.bonds(g.cells{i+1},1); % an array of all the vertices indices of the cell
vert = getRelativePosition(g,vidx,i); % the position of the vertices
nb=length(g.cells{i+1});
for j=1:nb
    next = mod(j,nb)+1;   % the next vertex
    n = norm(vert(j,:)-vert(next,:));
    adj_cell = g.bonds( g.cells{i+1}(j), 4 );
    if adj_cell ~= 0,
        ang_diff = g.polarity( adj_cell ) - g.polarity(i);
    else
        ang_diff = 0;
    end
    dp = dp + n*sign(cos(ang_diff))*sin(ang_diff);  
end
