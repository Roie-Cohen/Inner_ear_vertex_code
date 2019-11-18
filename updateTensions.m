function g = updateTensions(g)
% update the bond tension for each bond
% populations: 1 - general, 2 - SC, 3 - HC, 4 - Pillar, 5 - top border cell

for i=1:length(g.bonds(:,1))
    if g.bonds(i,1) ~= 0
        p = g.populations(g.bonds(i,3:4)); % the two types of seperated cells
        if sum(ismember([4 1], p)) == 2 || sum(ismember([4 2], p)) == 2 || sum(ismember([4 3], p)) == 2
            g.bonds(i,5) = 10; %10
        else
            if sum(ismember([5 2], p)) == 2 || sum(ismember([5 3], p)) == 2
                g.bonds(i,5) = 10;
            else
                if sum(p==[3 3])==2
                    g.bonds(i,5) = 5; 
                else
                    if sum(p==[2 2])==2
                        g.bonds(i,5) = 3; % 3
                    else
                        if sum(ismember([2 3], p)) == 2
                            g.bonds(i,5) = 1;
                        else
                            g.bonds(i, 5) = 1;
                        end
                    end
                end
            end
        end
    end
end

