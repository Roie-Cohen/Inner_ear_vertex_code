function g = differentiate(g, region)

    cc = 1:length(g.cells)-1;
    
    % seperating into two main populations
    g.populations = zeros(length(g.cells)-1,1);
    for i=cc
        c = cellCenter(g,i); 
        if abs(c(2)) < region
            g.populations(i) = 2; % the cells in the center
        else
            g.populations(i) = 1; % Hansen cells
%             g.areas(i) = 1.5*g.areas(i); % pillar cells are bigger
        end
    end
    
    % separating the hair cells and suporting cells
    perims = zeros(length(cc),1); % current perimeters of the cells
    for i=cc
        perims(i) = cellPerimeter(g, i);
    end
    for i = cc
        % lateral inhibition starting pattern
%         spec_cells = [97 99 88 102 79 104 50 26 52 29 54 43 74 107 78 93 57 71 46 48];
    % inner HCs pattern with pillar cells
       innerHCs = [25 26 16 31 45 47];
       spec_HCs = [innerHCs   61 51 53 66 57 71 82 77 75 85 80 108 102 104 106 109 115 100];
       spec_PCs = [37 38 39 40 29 41 42 43 56 59 60];
%        % for disord_lat14x14#1
%        innerHCs = [44 46 48 50 52 54 56];
%        spec_HCs = [innerHCs   85 73 89 91 78 107 95 110 124 126 136 134 120 132 117 115 142 162 145 152 163];
%        spec_PCs = [57:60 62:64 66:70];
       if (g.populations(i) == 2) && sum(spec_HCs == i) == 1
            g.populations(i) = 3; % hair cells
       end
       if (g.populations(i) == 2) && sum(spec_PCs == i) == 1
            g.populations(i) = 4; 
       end
    end
    
%     % convert Hensen cells to Pillar cells.
%     Hensen = find(g.populations == 1); % Hensen cells 
%     for i=1:length(Hensen)
%         c = Hensen(i);
%         b = g.cells{c+1}; % current cell bonds
%         ne = sum(g.bonds(b,3:4),2) - c; % current cell neighbors
%         ne_pop = g.populations(ne); % populations of neighbors
%         if ~isempty(find(ne_pop == 3,1)) || ~isempty(find(ne_pop == 2,1))
%             g.populations(c) = 4;
%             g.areas(c) = g.areas(c)*0.5;
%         end        
%     end 
    
    sym = 1;
%     if sym
        % 12x12
%         c = [3:3:144];
%         tmp = 4*[2:2:12];
%         tmp = [tmp-3 tmp-2 tmp-1 tmp];
%         c(tmp) = c(tmp) - 2;
%         % 6x6
%         c = [3:3:36];
%         tmp = 2*[2:2:6];
%         tmp = [tmp-1 tmp];
%         c(tmp) = c(tmp) - 2;

%     %12x12 spacious distribution
%     c = [14:2:24 37:2:47 62:2:72 85:2:95 110:2:120 133:2:143];
%         for i = cc
%             if ismember(i,c)
%                 g.populations(i) = 3;
%             end
%         end
%     end


%     % randomly distributing the hair cells and supporting cells
%     for i=cc
%         if g.populations(i)==2 && rand()<0.4
%             g.populations(i)=3;
%         end
%     end

end