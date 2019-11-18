nlat = 20;
rows = 12;
cols = 12;
for i=1:nlat
    g = disorderedLattice(rows, cols, 1);
    save(['lattices/disord_lat', num2str(rows), 'x', num2str(cols), '#', num2str(i)], 'g');
end