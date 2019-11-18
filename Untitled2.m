path = 'D:\Dima\project\sumer 2015\corent project\database\';
fileName = 'Series030_base_e15.5_12.4.mat';
filePath = strcat(path, fileName);

skel = a.SK;
[height, width] = size(skel);

% finding the vertices
verts = zeros(height, width);
nb = 0; % number of bonds
nv = 0; % number of vertices
for i=2:(height-1)
    for j=2:(width-1)
        if skel(i,j) == 1
            pix = [ [skel(i-1,j-1), skel(i-1,j), skel(i-1,j+1)];
                    [skel(i,j-1)  , skel(i,j)  , skel(i,j+1)  ];
                    [skel(i+1,j-1), skel(i+1,j), skel(i+1,j+1)] ];
            pixSum = sum(sum(pix));
            if pixSum >= 4
                verts(i,j) = 1;
                nb = nb + pixSum-1;
                nv = nv + 1;
            end
        end
    end
end

vidx = zeros(nv, 2);
[vidx(:,1), vidx(:,2)] = find(verts); % vertices indices

% finding the bonds
bonds = zeros(nb, 2);
b_counter = 0;
for i=1:nv
    rv = vidx(i,1);
    cv = vidx(i,2);
    pix = [ [skel(rv-1,cv-1), skel(rv-1,cv), skel(rv-1,cv+1)];
            [skel(rv  ,cv-1), 0            , skel(rv  ,cv+1)];
            [skel(rv+1,cv-1), skel(rv+1,cv), skel(rv+1,cv+1)] ];
    b = numel(find(pix));
    branches = zeros(b,2);
    [branches(:,1), branches(:,2)] = find(pix);
    for j=1:b
        b_counter = b_counter + 1;
        con_vert = findConnectingVert(skel, vidx, i, branches(j,:));
        bonds(b_counter, :) = [i, con_vert]; 
    end
end

figure;
hold on;
imshowpair(a.Im ,skel);
for i=1:nb
   if bonds(i,2)~=0
       vind1 = bonds(i,1);
       vind2 = bonds(i,2);
       r1 = vidx(vind1,:);
       r2 = vidx(vind2,:);
       line([r1(2),r2(2)],[r1(1),r2(1)]);
   end
end
% scatter(vidx(:,2),vidx(:,1), 5,'filled', 'r')



