function g = redistributeAreas(g)
% redistributes the preferable areas of the cells A0 (for example, after T2
% transition).

% a periodic lattice is scaled to 2pi on 2pi.
if g.bc==1
    Atot = (2*pi)^2;
else
    Atot = 363;
end
% number of cells of each type
Nh = sum(g.populations == 1 & ~g.dead);
Ns = sum(g.populations == 2 & ~g.dead);
NH = sum(g.populations == 3 & ~g.dead);
Np = sum(g.populations == 4 & ~g.dead);
Nt = sum(g.populations == 5 & ~g.dead);

% The areas are distributed so that the total area remains Atot
denom = Nh+g.fa(2)*Ns+g.fa(3)*NH+g.fa(4)*Np+g.fa(5)*Nt;
Ah = (g.fa(1)/(denom))*Atot;
As = (g.fa(2)/(denom))*Atot;
AH = (g.fa(3)/(denom))*Atot;
Ap = (g.fa(4)/(denom))*Atot;
At = (g.fa(5)/(denom))*Atot;

g.areas(g.populations == 1 & ~g.dead) = Ah;
g.areas(g.populations == 2 & ~g.dead) = As;
g.areas(g.populations == 3 & ~g.dead) = AH;
g.areas(g.populations == 4 & ~g.dead) = Ap;
g.areas(g.populations == 5 & ~g.dead) = At;

% Inner HCs (code is not the good...)
ypc = mean(g.verts(g.bonds(g.cells{find(g.populations==4,1)+1},1),2));
for c=1:length(g.cells)-1
    if g.populations(c) == 3 && mean(g.verts(g.bonds(g.cells{c+1},1),2)) < ypc
        g.areas(c) = As;
    end
end

end