function r_cm = cellCOM(g, i)
% calclates the center of mass of the cell
% i can inlude more than one index
% using a known formula for calculation of COM of a polygon
% area: A = 0.5*SUM(i=0:n-1)(x[i]*y[i+1]-x[i+1]*y[i]). where: x[n]==x[0]
if g.dead(i), return; end
r_cm = zeros(length(i) ,2);
for j=1:length(i)
    vidx=g.bonds(g.cells{i(j)+1},1); % an array of the vertices indices of the cell
    vert = getRelativePosition(g,vidx,i(j)); % the position of the vertices
    x = vert(:,1); xs = circshift(x,1);
    y = vert(:,2); ys = circshift(y,1);
    A = 0.5*( dot(x,ys) - dot(xs,y) );
    xcm = dot( x+xs , x.*ys - xs.*y )/(6*A);
    ycm = dot( y+ys , x.*ys - xs.*y )/(6*A);
    r_cm(j,:) = [xcm, ycm];
end

end