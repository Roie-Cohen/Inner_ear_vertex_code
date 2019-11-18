function Rcm = COM(g, pop)
% calculates the center of mass of cell type 'pop'
% the average weight for each cell is relative to cell size
c = find(g.populations == pop);
nc = length(c);
tot = [0 0];
tot_area = 0;
for i=1:nc
    cent = cellCenter(g,c(i));
    A = cellarea(g,c(i));
    tot_area = tot_area + A;
    tot = tot + A*cent;
end
Xcm = tot(1)/tot_area;
Ycm = tot(2)/tot_area;
Rcm = [Xcm Ycm];

end
