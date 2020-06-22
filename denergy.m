function dE=denergy(v,gin)

dE=zeros(2*length(gin.verts),1);
g=insertverts(v,gin);
for i=1:length(g.cells)-1
    if ~g.dead(i)
        
        dEarea = dHarea(g,i);
        dElength = dHlength(g,i);
        dEsteric = dHrepulsion(g, i);
        dEshear = dHshear(g, i, 2.4);
        dEcomp = dHcompression(g, i);
        
        dEi = dEarea + dElength + dEsteric + dEcomp + dEshear;
        dE = dE+dEi;
    end
end

dE=dE';
end