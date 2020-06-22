function g = updateParameters(g)
% updates cells' mechanical parameters according to cell type
% populations: 1 - general, 2 - SC, 3 - HC, 4 - Pillar, 5 - top border cell

if isfield(g, 'populations')
    nc = length(g.cells) - 1;
    % compressibility and contractility
    for i=1:nc
        g.cmpr(i) = g.alpha(g.populations(i));
        g.ctrc(i) = g.Gamma(g.populations(i));
    end
    
    % update the bond tension for each bond
    for i=1:length(g.bonds(:,1))
        if g.bonds(i,1) ~= 0
            p = g.populations(g.bonds(i,3:4)); % the two types of seperated cells
            g.bonds(i,5) = g.tension(p(1),p(2));
        end
    end

    
end
end