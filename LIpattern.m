function g = LIpattern(g, yreg)
% get's disordered non-differentiated lattice and creates a pattern similar to 
% lateral inhibition pattern.
% yreg is the region to be differentiated, from yreg(1) to yreg(2)

nc = length(g.cells)-1;
cperm = randperm(nc); % random permutation of the cells
g.populations = 2*ones(nc, 1); % initiating every cell to start as SC
for i=1:nc
    c = cperm(i); % current cell index
    if g.dead(c) % if the cell is dead then the loop skips this cell
        continue;
    end
    cpos = cellCenter(g, c);
    if cpos(2)<yreg(1) || cpos(2)>yreg(2)   % if the cell isn't in the right region then the loops skips it
        continue;
    end
    cneigh = g.bonds(g.cells{c+1},4); % nieghboring cells to 'c'
    npop = g.populations(cneigh);   % neighbors population
    if isempty(find(npop==3,1))
        g.populations(c) = 3;
    end
end % end of for loop

% end of function
end
    
    