function frac = surroundingHCsFrac(g,i)
if ~g.dead(i)
    perim=0;
    HCs_l = 0;
    vidx=g.bonds(g.cells{i+1},1);
    vert = getRelativePosition(g,vidx,i);
    nb=length(vidx);
    for j=1:nb
        n = norm(vert(j,:)-vert(mod(j,nb)+1,:));
        perim=perim+n;
        HCs_l = HCs_l + n*(g.populations(g.bonds(g.cells{i+1}(j),4))==3);
    end
    frac = HCs_l/perim;
end
end