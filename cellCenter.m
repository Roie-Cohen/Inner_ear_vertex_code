function x = cellCenter(g, i)
% calclates the center of mass of the cell
for j=1:length(i)
    vidx=g.bonds(g.cells{i(j)+1},1); % an array of the vertices indices of the cell
    vert = getRelativePosition(g,vidx,i(j)); % the position of the vertices
    p = polyshape(vert(:,1),vert(:,2));
    [xcm, ycm] = centroid(p);
    x(j,:) = [xcm, ycm];
%     x(j,:) = mean(vert);
end

end