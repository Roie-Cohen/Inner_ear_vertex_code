function A = average_area(g, cells)
% cells: cells indices
    cells = cells(~g.dead(cells));
    nc = length(cells);
    areas = zeros(nc, 1);
    for i=1:nc
        areas(i) = cellarea(g,cells(i));
    end
    A = mean(areas);
end