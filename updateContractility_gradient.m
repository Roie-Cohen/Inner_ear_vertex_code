function g = updateContractility(g)
% updates the contactility for each cell
% 18/07/2018:
% contractility of is set higher closer to the pillar cells
% contractility is reduced with an exponential gradient

sigma = 2; % exponential decay length

nc = length(g.cells)-1;

PCs = find(g.populations == 4);
PCs_pos = cellCenter(g, PCs);
Lin_fit = polyfit(PCs_pos(:,1), PCs_pos(:,2), 1); % line defining the pillar cells axis
m = Lin_fit(1); %mx+n
n = Lin_fit(2);

for i=1:nc
    p = g.populations(i);
    if p == 2 || p == 3
        cpos = cellCenter(g, i);
        x = cpos(1); y = cpos(2);
        D = (y-m*x-n)/sqrt(m^2+1); % distance from the pillar line
        g.ctrc(i) = 1+2*exp(-D/sigma)*(D>0);
    else
        g.ctrc(i) = 1;
    end
end

end