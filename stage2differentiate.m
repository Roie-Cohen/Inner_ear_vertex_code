function g = stage2differentiate(g)
% gets a g structure of a lattice after stage1
% adds the top boundary cells and seperates SCs and general cells.

cc = 1:length(g.cells)-1; % total number of cells
% changes every SC (population=2) to general cell it doesn't touch a HC
for i = cc
    if (g.populations(i) == 2) && isempty(find(g.populations(g.bonds(g.cells{i+1},4))==3,1))
        g.populations(i) = 1;
    end
end

for i = cc
    if ~g.dead(i)
        % changes every SC in the upper half of the lattice to top-border
        % cell if it borders a general cell.
        cpos = cellCOM(g,i);
        if cpos(2)>0 && (g.populations(i) == 2) &&  ~isempty(find(g.populations(g.bonds(g.cells{i+1},4))==1,1))
            g.populations(i) = 5;
        end
        % changes every SC in the lower half of the lattice to temporary
        % type cell if it borders a general cell.
        if cpos(2)<0 && (g.populations(i) == 2) && ~isempty(find(g.populations(g.bonds(g.cells{i+1},4))==1,1))
            g.populations(i) = 11; % temp type
        end
    end
end
% changes all temporary types to general cells
g.populations(g.populations == 11) = 1;

for i = cc
    if ~g.dead(i)
        % changes a general cell that touches 2 HCs but no SC to SC
        if (g.populations(i) == 1) &&  isempty(find(g.populations(g.bonds(g.cells{i+1},4))==2,1)) && sum(g.populations(g.bonds(g.cells{i+1},4))==3) ==2
            g.populations(i) = 2;
        end
        % changes a top border cell to SC if it touches 3 HCs or more
        if (g.populations(i) == 5) && sum(g.populations(g.bonds(g.cells{i+1},4))==3) >= 3
            g.populations(i) = 2; 
        end
    end
end

for i = cc
    if ~g.dead(i)
        % changes a general cell in the upper half of the lattice to top border cell if it touches a SC
        cpos = cellCOM(g,i);
        if cpos(2)>0 && (g.populations(i) == 1) &&  ~isempty(find(g.populations(g.bonds(g.cells{i+1},4))==2,1))
            g.populations(i) = 5;
        end
    end
end


end