% get several statistics about the lattice

neigh = zeros(length(g.cells)-1,3);
outer = zeros(length(g.cells)-1,1);
for i=1:length(g.cells)-1
    for j=1:length(g.cells{i+1})
        border_cells = g.bonds(g.cells{i+1}(j), 3:4);
        if ismember(0, border_cells)
            outer(i) = 1;
            break;
        end
        np = g.populations(border_cells(find(g.populations(i) ~= g.populations(border_cells),1)));
        if g.populations(border_cells(1))==g.populations(border_cells(2))
            np = g.populations(border_cells(1));
        end
        neigh(i, np) = neigh(i, np) + 1; 
    end
end

HCs = sum(g.populations(~outer) == 3);
SCs = sum(g.populations(~outer) == 2);
others = sum(g.populations(~outer) == 1);

disp('**************************************');
disp(strcat('**   Hair cells:', 9,9,9 ,num2str(HCs),9,9, '**'));
disp(strcat('**   Supporting cells:', 9,9 ,num2str(SCs),9,9, '**'));
disp(strcat('**   Other types of cells:', 9,num2str(others),9,9, '**'));
disp('**************************************');

HSr = HCs/SCs;
disp(strcat('**   HCs to SCs ratio:', 9, 9,num2str(round(HSr,2)),9, '**'));
