function dE=denergy(v,gin)

dE=zeros(2*length(gin.verts),1);
g=insertverts(v,gin);
for i=1:length(g.cells)-1
    if(g.dead(i) ==0)
        if isfield(g,'populations') && g.populations(1)~=0
            vdv = dHvdv2(g, i);
            pcs = dHpcs3(g, i);
            flow = dHflow4(g, i, 2.4);
        else
            vdv = 0;
            pcs = 0;
            flow = 0;
        end
     dEi = dHarea(g,i) +dHlength(g,i) + vdv + pcs + flow;
     
     dE = dE+dEi;

    end
end

dE=dE';
end