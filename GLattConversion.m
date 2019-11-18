 function g = GLattConversion(Cells,Vertices) % Converts Voronoi-structures into GLattice structures. 
 
    
    NCells = length(Cells);
    NVertices = length(Vertices);
    g = struct('cells',[],'bonds',[],'verts',[]);
    
    
    % Populate Vertices
    g(1).verts = zeros(NVertices,3);
    g(1).verts(:,1:2) = Vertices;
    
    
    % Make sure cells are sorted! 
    for I = 1 : NCells
        [ClockSorted,Index] = ClockWiseSort(Cells{I}(1:(length(Cells{I})-1)),Vertices);
        Cells{I} = [Cells{I}(Index),Cells{I}(Index(1))];
    end
    
    % Convert cell entries into bonds and note: 
    % convention for bond=[v1 v2 c1 c2], looking from v1 to v2 -> c1= left cell
    % and c2->right cell
    [g(1).cells,g(1).bonds] = BondsConversion(Cells,NVertices);
    
    
    
    function [cells,bonds] = BondsConversion(Cells,NVertices)
        % determine NBonds 
        NBonds = 0;
        for i = 1 : NCells
            NBonds = NBonds + length(Cells{i}) - 1;
        end

		%% bonds = zeros(NBonds,4);
%% added
        bonds = zeros(NBonds,5);
        Counter = 1;
        cells = cell(1,NCells+1);%NCells+1
        cells{1} = [];%zeros(1,size(ConvexHull,1)-1);% The outside cell, i.e. cell0 in c-code.
 
        S=sparse(NVertices,NVertices);
        for i = 1 : NCells
            for j = 1 : (length(Cells{i})-1)             
                v1 = Cells{i}(j);
                v2 = Cells{i}(j+1);
                S(v1,v2) = i;
            end
        end       
        
        for i = 1 : NCells
            
            cells{i+1} = zeros(1,length(Cells{i})-1);
            for j = 1 : (length(Cells{i})-1)
                v1 = Cells{i}(j);
                v2 = Cells{i}(j+1);
                % find c1 and c2:
                
                if S(v2,v1) ~= 0
                    % inside bond
                    cells{i+1}(j) = Counter; %i+1 
                    bonds(Counter,1:4) = [Cells{i}(j) Cells{i}(j+1) i S(v2,v1)];
                   
                else
                    % outside bond
                    bonds(Counter,1:4) = [Cells{i}(j) Cells{i}(j+1) i 0];
                    cells{i+1}(j) = Counter; %i+1
                    Counter = Counter + 1;
                    bonds(Counter,1:4) = [Cells{i}(j+1) Cells{i}(j) 0 i];
                    cells{1} = [cells{1},Counter];
                end
                Counter = Counter + 1;
            end
        end
       
        % Finally cells{1} needs some resorting: 
        ToSort = bonds(cells{1},1:2);
        NewCell = zeros(size(cells{1}));
        NewCell(1) = 1;
        for i = 1 : (length(cells{1})-1),
            NewCell(i+1) = find(ToSort(:,2) == ToSort(NewCell(i),1));
        end
        NewCell = NewCell(length(NewCell):-1:1);
        cells{1} = cells{1}(NewCell);
        
    end

    function [ClockSorted,Index] = ClockWiseSort(NewCell,NewVertices)
        NewCellLength = length(NewCell);
        %   first we need to find the center of mass of the vertices
        if NewCellLength == 1,
            Center = NewVertices(NewCell,:);
        else
            Center = sum(NewVertices(NewCell,:))/NewCellLength;
        end
        %   Now take the NewVertices, and subtract the Center from each entry
        MovedVertices = NewVertices(NewCell,:);
        MovedVertices(:,2) = MovedVertices(:,2) - Center(2);
        MovedVertices(:,1) = MovedVertices(:,1) - Center(1);
        [Sorted,Index] = sort(cart2pol(MovedVertices(:,1),MovedVertices(:,2))); 
        % Now we actually sort counter clock, just invert.
        Index = Index(length(Index):-1:1);
        ClockSorted = NewVertices(NewCell(Index),:); 
    end

%%set the tension parameter
g.bonds(:,5)=0.1;
    
end