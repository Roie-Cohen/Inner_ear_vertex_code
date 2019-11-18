figure(11)
hold on
for c=1:length(g.cells)-1
    vidx = g.bonds(g.cells{c+1}, 1);
    vorg = g.verts(vidx, 1:2);
    vrel = getRelativePosition(g, vidx, c);
    dv = vorg-vrel;
    ri = find(sum(abs(dv),2) ~= 0);
    reloc = vidx(ri);
    for i=1:length(reloc)
%         ind = reloc(i);
        vpos = vrel(ri(i),:);
        nind = length(g.verts(:,1))+1; % new vertex index
        g.verts = [g.verts; [vpos 0]];
        % update bonds
        g.bonds(g.cells{c+1}(ri(i)), 1) = nind;
        g.bonds(g.cells{c+1}(ri(i)), 4) = 0;
        g.bonds(g.cells{c+1}(mod(ri(i)-2,length(vidx))+1), 2) = nind;
        g.bonds(g.cells{c+1}(mod(ri(i)-2,length(vidx))+1), 4) = 0;
    end
    scatter(g.verts(reloc, 1), g.verts(reloc, 2));
end