function g = getGfromCIA(rotate_rad)
% selecting and loading the cia (cells image analysis) 
db_path = 'Q:\Users\Roie\cia_example\';
% filename = 'base 2_e17.5_28.1.mat';
filename = 'base 3_e17.5_28.1.mat';
filename = uigetfile(db_path);
load([db_path filename]);


% remove bonds (sides) with no vertices (errors in segmentation)
count = 0;
for i=1:length(a.cia.sides)
    if isempty(a.cia.sides(i-count).vert)
        a.cia.sides(i-count) = [];
        for j=1:length(a.cia.cells) % update the indices of the sides in the cells structure
            a.cia.cells(j).sides(a.cia.cells(j).sides >= i-count) = a.cia.cells(j).sides(a.cia.cells(j).sides >= i-count) - 1;
        end
        count = count + 1;
    end
end


% start to build the 'g' structure
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
g.dead =zeros(nc,1);
g.cells = cell(1,nc+1);
for i=1:nc+1
    g.cells{i} = a.cia.cells(i).sides;
    if isempty(g.cells{i})
        g.dead(i-1) = 1;
    end
end


% corrections for single vertex bonds
svb = find( g.bonds(:,1)-g.bonds(:,2) == 0 ); % single vertex bonds
for i=1:length(svb)/2
    b1 = svb(i);
    for k=i+1:length(svb)
        if sum(ismember(g.bonds(b1,3:4),g.bonds(svb(k),3:4)))==2 %the two bonds seperate the same cells
            if k~= i+1
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
%     for j=1:2
%         g.cells{c(j)}(g.cells{c(j)} == b1+i-1) = b1;
%         g.cells{c(j)}(g.cells{c(j)} == b2+i-1) = [];
%     end
    for j=1:length(g.cells)
        g.cells{j}(g.cells{j}==b2) = b1;
        g.cells{j} = unique(g.cells{j});
        g.cells{j}(g.cells{j}>b2) = g.cells{j}(g.cells{j}>b2) - 1;
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

global g1
g1 = g;

% correction of adjecent vertices (vertices with one pixel separation).
% deals with groups of doubles, triplets and quadruplets vertices.
% replacing each group of close vertices with a single vertex.
v = [];
for i=1:length(g.verts)
    for j=1:length(g.verts)
        if norm(g.verts(i,1:2) - g.verts(j,1:2)) < 2 && i~=j
            v = [v; [i j]];
        end
    end
end
% disp(v)
tot_verts = [];
while ~isempty(v)
    group_start = v(1,:);
    group_ind = find(v(:,1) == group_start(1) | v(:,2) == group_start(1) | v(:,1) == group_start(2) | v(:,2) == group_start(2));
    vgroup = v(group_ind,:);
    vgroup = unique(vgroup);
    tot_verts = [tot_verts; vgroup];
    % convert the group to a single vertex
    new_vind = length(g.verts(:,1)) + 1;
    new_pos = sum(g.verts(vgroup,:))/length(vgroup);
    g.verts(new_vind,:) = new_pos;  % adding the new vertex
    % linking the relavent bonds to the new vertex
    for i=1:length(vgroup)
        for j=1:2
            g.bonds(g.bonds(:,j) == vgroup(i), j) = new_vind;
        end
    end
    % remove this group from v
    for i=1:length(vgroup)
        v(v(:,1)==vgroup(i) | v(:,2)==vgroup(i),:) = [];
    end
end
% remove old vertices from g and redefine the vertices in g.bonds
tot_verts = sort(tot_verts,'descend');
for i=1:length(tot_verts)
    g.verts(tot_verts(i),:) = [];
    for j=1:2
        g.bonds(g.bonds(:,j) > tot_verts(i), j) = g.bonds(g.bonds(:,j) > tot_verts(i), j) - 1;
    end
end
% remove single vertex bonds created by combining the vertices groups
svb = find(g.bonds(:,1)-g.bonds(:,2)==0);
svb = sort(svb,'descend');
for i=1:numel(svb)
    bcells = g.bonds(svb(i),3:4);
    for j=1:2
        if sum(g.bonds(:,3) == bcells(j) | g.bonds(:,4) == bcells(j)) == 1
            g.dead(bcells(j)) = 1;
        end
        g.cells{bcells(j)+1}(g.cells{bcells(j)+1}==svb(i)) = [];
    end
    g.bonds(svb(i),:)=[];
    for j=1:length(g.cells)
        g.cells{j}(g.cells{j}>svb(i)) = g.cells{j}(g.cells{j}>svb(i))-1;
    end
end

        
global g2
g2 = g;

% organize the bonds of each cell in a clockwise direction
for i=2:nc+1
    if ~g.dead(i-1)
        bb = g.cells{i};
        bb_new = bb;
%         disp(i)
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

% rotate the choclea in 'rotate_rad' radians
if nargin == 1
    g = rotateLattice(g, rotate_rad, 0);
end

% % specifically for 'base 2_e17.5_28.1.mat'
% others = [1:12 54 56 62 63 64 66:71];
% g.populations(others) = 1;

% specifically for 'e16.5_samp1apex-a.mat'
others = [2 31 40 55 66 81 98 120 143 138 8 17 29 45 58 72 84 97 110 126 139 4 20 28 44 61 75 89 101 113 121 133 145 14 21 38 52 65 80 92 3 11 22 39 47 62 67 74 85 100 112 108 118 127 134 144 25 46 56 76 86 99 106 117 131 140 152 18 34 36 63 79 95 115 130 149];
g.populations(others) = 1;

end