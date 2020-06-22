function [ng  nkilled]= killSmallCells(g,thresh,crit)
nkilled = 0;
if (nargin<3)
    crit = 1;
end
ng  =g;
nbcell = length(ng.cells)-1;
cand = find(ng.dead==0);
for ci=1:length(cand)
    c = cand(ci);
    if (cellarea(ng,c)<thresh*ng.areas(c))
        while (length(ng.cells{c+1})>3)
            if(crit==1) %% sort by boundary length
                vidx=g.bonds(ng.cells{c+1},1);
                vert = getRelativePosition(ng,vidx,c);
                nb=length(vidx);
                len = zeros(nb,1);
                for j=1:nb
                    len(j) = norm(vert(j,1:2)-vert(mod(j,nb)+1,1:2));
                end
                [ord, id] =sort(len);
            elseif (crit==2), %% sort by neighbouring notch level
                neib =  g.bonds(ng.cells{c+1},4);
%                 [ord, id] =sort(g.level(neib,2));
                [ord, id] =sort(g.delta(neib));
            elseif (crit==3); %% sort by orientation
                vidx=g.bonds(ng.cells{c+1},1);
                vert = getRelativePosition(ng,vidx,c);
                nb=length(vidx);
                ori = zeros(nb,1);
                for j=1:nb
                    ori(j) = abs(cos(atan2(vert(j,2)-vert(mod(j,nb)+1,2),vert(j,1)-vert(mod(j,nb)+1,1) )));
                end
                [ord, id] =sort(ori);
            end
            ng.transitionedBonds(ng.transitionedBonds(:,1) == ng.cells{c+1}(id(1)),:) = [];
            ng = T1transition(ng,ng.cells{c+1}(id(1)));
        end
        ng = T2transition(ng,c);
        nkilled = nkilled+1;
    end
end
