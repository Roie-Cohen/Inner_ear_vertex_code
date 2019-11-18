function g=findTransitions(g,eps,probT1,probT2)
%% generates T1 transition with probability  probT1 if edge is  shorter than eps
%% if cell has less than 4  edges,  generates a T2 transition is probability probT2
cand = randperm(size(g.bonds,1)); %% random  selection order
while ~isempty(cand)
    bo = cand(1);
    if (g.bonds(bo,1)==0)
        cand = cand(2:end);
        continue;
    end
    if (g.bc==2)
        if (~isempty(find(g.fixed_verts(:,1) == g.bonds(bo,1), 1)) || ~isempty(find(g.fixed_verts(:,1) == g.bonds(bo,2), 1)) ),
            cand = cand(2:end);
            continue;
        end
    end
    
    v1 = g.verts(g.bonds(bo,1),:);
    v2 = g.verts(g.bonds(bo,2),:);
    vec = v1-v2;
    
    if(length(g.cells{g.bonds(bo,3)+1})>3 && length(g.cells{g.bonds(bo,4)+1})>3)
        len = norm(vec);
        if(len<eps)
            if isfield(g,'populations') && g.populations(1)~=0
                edge_cells = find_bond_edge_cells(g, bo);
                edge_pop = g.populations(edge_cells);
                % if the bond seperates HC:HC there's no transition
                if  sum(edge_pop==[3;3]) ~= 2
                    if rand() < probT1
                        g = T1transition(g,bo);
                    end
                end

            else
                if rand()<probT1
                    g = T1transition(g,bo);
                end
            end
        end
        
    else
        c1 = g.bonds(bo,3);
        c2 = g.bonds(bo,4);
        torem = [];
        if(length(g.cells{c1+1})==2) % 3
            if(rand()<probT2),
                [g torem]=T2transition(g,c1);
                g = redistributeAreas(g);
            end
        end
        if(length(g.cells{c2+1})==2) % 3
            if(rand()<probT2),
                [g torem]=T2transition(g,c2);
                g = redistributeAreas(g);
            end
        end
    end
    cand = cand(2:end);
end

