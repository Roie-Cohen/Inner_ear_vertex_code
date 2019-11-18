function g = updateCellsCenter(g)
% updates the centroid position of the cells

for i=1:length(g.cells)-1
    if ~g.dead(i) && (g.populations(i)==3 || g.populations(i)==4)
        g.centroid(i,:) = cellCOM(g,i);
    end
end

end