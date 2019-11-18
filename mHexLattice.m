function [HexCells,Vertices] = mHexLattice(NRows,NCols,rnd)

    if nargin < 1
        NRows = 5;  % NRows should be odd 
        NCols = 10;
    end
    NCells = NRows*NCols;
    CellsPositions = zeros(NCells,2);
    Neighbours = cell(NCells,1);
    for i = 1 : NCells
        Neighbours{i} = [i-1,i-NCols,i-NCols-1,i+NCols,i+NCols+1,i+1];
    end
    for i = 1 : NRows
        for j = 1 : NCols
            CellsPositions((i-1)*NCols+j,2) = i*3/2 +2*rnd*rand()-rnd;
            CellsPositions((i-1)*NCols+j,1) = j*sqrt(3)+2*rnd*rand()-rnd;
        end
    end
    for i = 1 : 2 : NRows
        CellsPositions(((i-1)*NCols+1):i*NCols,1) = CellsPositions(((i-1)*NCols+1):i*NCols,1) + sqrt(3)/2;
    end
    
    % Voronoi Diagramme
    [Vertices,Cells] = voronoin(CellsPositions);
    CellLengths = zeros(length(Cells),1);
    for i = 1 : length(Cells)
        CellLengths(i) = length(Cells{i});
    end
    Hexagons = find(CellLengths == 6); % This is important to adress the neighbours again!
    Correspondence = zeros(length(Cells),1);
    for i = 1 : length(Cells)
        if ~isempty(find(Hexagons == i))
            Correspondence(i) = find(Hexagons == i);
        end
    end
    HexCells = cell(length(Hexagons),1);
    for i = 1 : length(HexCells)
        HexCells{i} = Cells{Hexagons(i)};
    end
     %Vertices(Vertices == inf) = 0;
    if(0),
    HexNeighbours = cell(length(HexCells),1);
    NCols = NCols - 2
    NRows = NRows - 2
    for i = 1 : NRows
        for j = 1 : NCols
            N = NCols*(i-1) + j;
            if mod(i,2) == 0 % In this case i is even
                if j ~= 1 && j ~= NCols
                    HexNeighbours{N} = [N-1,N+1,N+NCols,N+NCols+1,N-NCols,N-NCols-1]; 
                elseif j == 1
                    HexNeighbours{N} = [N+NCol-1,N+1,N+NCols,N+NCols+1,N-NCols,N-NCols-1];
                elseif j == NCols
                    HexNeighbours{N} = [N-1,N-NCols+1,N+NCols,N+NCols+1,N-NCols,N-NCols-1];
                end
            else
                disp('need even number of rows');
            end
        end
    end
    NCells = length(HexCells);
    for i = 1 : length(HexNeighbours)
        HexNeighbours{i} = HexNeighbours{i}(find(HexNeighbours{i} > 0));
        HexNeighbours{i} = HexNeighbours{i}(find(HexNeighbours{i} <= NRows*NCols));
        if length(HexNeighbours{i}) < 6
            temp = (NCells+1)*ones(1,6);
            temp(1:length(HexNeighbours{i})) = HexNeighbours{i};
            HexNeighbours{i} = temp;
        end
    end
    end
end