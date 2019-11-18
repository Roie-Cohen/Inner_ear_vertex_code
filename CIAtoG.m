
% selecting and loading the cia (cells image analysis) 
db_path = 'Q:\Users\Roie\database\';
filename = 'base 2_e17.5_28.1.mat';
% filename = 'base 3_e17.5_28.1.mat';
% filename = uigetfile(db_path);
load([db_path filename]);

g = struct('cells',[],'bonds',[],'verts',[]);

nv = size(a.cia.vertices, 2);
g.verts = zeros(nv, 3);
for i=1:nv
    g.verts(i,:) = [a.cia.vertices(i).coord 0];
end

nb = size(a.cia.sides, 2);
g.bonds = zeros(nb, 5);
for i=1:nb
    g.bonds(i,1:2) = a.cia.sides(i).vert;
    g.bonds(i,3:4) = a.cia.sides(i).cells - 1;
end

nc = size(a.cia.cells, 2) - 1;
g.cells = cell(1,nc+1);
for i=1:nc+1
    g.cells{i} = a.cia.cells(i).sides;
end

% corrections for single vertex bonds
svb = find( g.bonds(:,1)-g.bonds(:,2) == 0 ); % single vertex bonds
for i=1:length(svb)/2
    b1 = svb(i);
    for k=i+1:length(svb)
        if sum(ismember(g.bonds(b1,3:4),g.bonds(svb(k),3:4)))==2 %the two bonds seperate the same cells
            if k~= i+1
               tmp = g.bonds(k,:);
               g.bonds(k,:) = g.bonds(i+1,:);
               g.bonds(i+1,:) = tmp;
               tmp = svb(k);
               svb(k) = svb(i+1);
               svb(i+1) = tmp;
            end
            break;
        end
    end
    b2 = svb(i+1);
    g.bonds(b1,2) = g.bonds(b2,1);
    c = g.bonds(b2, 3:4)+1;
    g.bonds(b2,:) = [];
    for j=1:2
        g.cells{c(j)}(find(g.cells{c(j)} == b1+i-1)) = b1;
        g.cells{c(j)}(find(g.cells{c(j)} == b2+i-1)) = [];
    end
end
nb = nb - length(svb)/2;

%  corrections for missing short bonds
skip = 0;
for i=1:nv
    if numel(a.cia.vertices(i).sides) == 2
        if skip == 1
            skip = 0;
            continue;
        end
        nb = nb + 1;
        cc1 = a.cia.vertices(i).cells;
        cc2 = a.cia.vertices(i+1).cells;
        count = 0;
        c = [0 0];
        for j=1:length(cc1)
            if ~isempty(find(cc2 == cc1(j),1))
                count = count + 1;
                c(count) = cc1(j);
            end
        end
        g.bonds = [g.bonds; [i i+1 c(1)-1 c(2)-1 0]];
        g.cells{c(1)} = [g.cells{c(1)} nb];
        g.cells{c(2)} = [g.cells{c(2)} nb];
        skip = 1;
    end
end

% very specific corrections for base3_e17.5_28.1.mat
if strcmp(filename, 'base 3_e17.5_28.1.mat')
    g.cells{3}(2) = 3;
    g.cells{46} = [116,119,125,127,130,159,166,219];
    g.cells{63} = [164,168,173,174,200,219];
end

% % draw the cells one by one
% for i=1:length(g.cells)
%     for j=1:length(g.cells{i+1})
%         line(g.verts(g.bonds(g.cells{i+1}(j),1:2),1), g.verts(g.bonds(g.cells{i+1}(j),1:2),2));
% %         pause(0.5)
%     end
%     pause(0.1)
% end


% organize the bonds of each cell in a clockwise direction
for i=2:nc+1
    bb = g.cells{i};
    bb_new = bb;
    vnext = g.bonds(bb(1),1);
    for j=2:length(bb)
        bidx = [find(g.bonds(bb,1) == vnext) find(g.bonds(bb,2) == vnext)];
        bnext = bb(bidx);
        bnext = bnext(find(bnext~=bb_new(j-1),1));
        bb_new(j) = bnext;
        vn = g.bonds(bnext, 1:2);
        vnext = vn(find(vn ~= vnext,1));
    end
    % now to ensure clockwise and not anti-clockwise
    v = g.bonds(bb_new, 1);
    if ispolycw(g.verts(v,1), g.verts(v,2) ) == 0
        bb_new = wrev(bb_new); %reverse the vector
    end
    g.cells{i} = bb_new;
end

% refining the g.bonds array to be in the right format
% g.bonds(b,1) is the first vertex of the bond b in the clockwise direction
% g.bonds(b,3) is the cell for the clockwise direction was taken for
newb = size(g.bonds, 1) + 1;
for i=1:length(g.cells)-1
    b = g.cells{i+1};
    for j=1:length(b)
        b1 = b(j);
        b2 = b( mod(j,length(b))+1 );
        if ismember(g.bonds(b1,1), g.bonds(b2,1:2))
            g.bonds = [g.bonds; [g.bonds(b1,2), g.bonds(b1,1), g.bonds(b1, 3:5)] ];
            g.cells{i+1}(j) = newb;
            b1 = newb;
            newb = newb + 1;
        end
        if g.bonds(b1,3) ~= i
            g.bonds(b1, 3:4) = wrev(g.bonds(b1, 3:4));
        end
    end
end

% rescale the lattice to scale 10
rescale = 20;
xscale = max(g.verts(:,1)) - min(g.verts(:,1));
yscale = max(g.verts(:,2)) - min(g.verts(:,2));
g.verts(:,1) = g.verts(:,1)*rescale/xscale;
g.verts(:,2) = g.verts(:,2)*rescale/yscale;


g.xboundary = zeros(length(g.cells)-1,2);
g.yboundary = zeros(length(g.cells)-1,2);
g.bc = 0;
g.scale = eye(2);
g.dead =zeros(length(g.cells)-1,1);

% populations
g.populations = zeros(nc,1);
for i=1:nc
    switch a.cia.cells(i+1).type
        case 0
            g.populations(i) = 2;
        case 1
            g.populations(i) = 3;
    end
end

% areas. HCs are larger than SCs by a factor af.
tot_area = 0;
for i=1:length(g.cells)-1
    tot_area = tot_area + cellarea(g,i);
end
af = 1.5;
SCs = sum(g.populations == 2);
HCs = sum(g.populations == 3);
g.areas = ones(length(g.cells)-1,1)*tot_area/(HCs*af + SCs);
g.areas(g.populations == 3) = g.areas(g.populations == 3)*af;


% % specifically for 'base 2_e17.5_28.1.mat'
% others = [1:12 54 56 62 63 64 66:71];
% g.populations(others) = 1;

% specifically for 'base 3_e17.5_28.1.mat'
% others = [1:12 54 56 62 63 64 66:71];
% g.populations(others) = 1;

