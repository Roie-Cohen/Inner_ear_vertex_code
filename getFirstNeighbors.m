function fn = getFirstNeighbors(ipos, npos)
% Gets the position of the point of interest 'ipos' and the positions of all 
% neighboring points 'npos' in [x1 y1; x2 y2; ...] format.
% The function returns the indices of the first neighbors according to
% Voronoi tessalation. The returned indices are the indices of the points in 'npos'.

[~, cc] = voronoin([ipos; npos]);

icc = cc{1}; % the vertices of the point of interest
fn = zeros(length(icc)-1, 1); % initializing number of first neighbors
count = 1;
for c = 2:length(cc)
    if sum(ismember(cc{c}, icc)) ~= 0
        fn(count) = c - 1;
        count = count + 1;
    end
end

end