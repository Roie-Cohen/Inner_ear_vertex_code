function g = convert2pillar(g)
% Find all the SCs (2) that touches Hancen (1) and coverts them to pillar
% cells (4).
% If a pillar cell is pushed out (not touching any SC or HC) it is
% converted to Hancen (1)
% If a certain time has passed from the start of the simulation, a HC touching a Hancen will turn to pillar.
global timer

%     SCs = find(g.populations == 2); % SCs 
%     for i=1:length(SCs)
%         c = SCs(i);
%         b = g.cells{c+1}; % current cell bonds
%         ne = sum(g.bonds(b,3:4),2) - c; % current cell neighbors
%         ne_pop = g.populations(ne); % populations of neighbors
%         if ~isempty(find(ne_pop == 1,1))
%             g.populations(c) = 4;
%             g.areas(c) = g.areas(c)*0.5;
%         end
%     end
    
%     PCs = find(g.populations == 4); % Pillar cells 
%     for i=1:length(PCs)
%         c = PCs(i);
%         b = g.cells{c+1}; % current cell bonds
%         ne = sum(g.bonds(b,3:4),2) - c; % current cell neighbors
%         ne_pop = g.populations(ne); % populations of neighbors
%         if isempty(find(ne_pop == 2,1)) && isempty(find(ne_pop == 3,1))
%             g.populations(c) = 1;
%             g.areas(c) = g.areas(c)*2;
%         end        
%     end    
    
%     HCs = find(g.populations == 3); % HCs
%     for i=1:length(HCs)
%         c = HCs(i);
%         b = g.cells{c+1}; % current cell bonds
%         ne = sum(g.bonds(b,3:4),2) - c; % current cell neighbors
%         ne_pop = g.populations(ne); % populations of neighbors
%         if ~isempty(find(ne_pop == 4,1))
%             g.areas(c) = 2*average_area(g, find(~g.dead));
%         end        
%     end  
    
%     Hancen = find(g.populations == 1); % Hancen cells 
%     for i=1:length(Hancen)
%         c = Hancen(i);
%         b = g.cells{c+1}; % current cell bonds
%         ne = sum(g.bonds(b,3:4),2) - c; % current cell neighbors
%         ne_pop = g.populations(ne); % populations of neighbors
%         if timer>500 && ~isempty(find(ne_pop == 3,1))
%             g.populations(c) = 4;
%             g.areas(c) = g.areas(c)*0.5;
%         end        
%     end 
end