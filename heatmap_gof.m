function heatmap_gof(param_name)

switch param_name
    case 'area'
        ps = [2 4 6 9 13 18];
    case 'repulsion'
        ps = [0 0.2 0.5 1 1.5 3 5];
    case 'tension'
        ps = [0.05 0.1 0.2 0.3 0.45 0.6 0.8];
    case 'shear'
        ps = [0 0.1 0.2 0.4 0.6 0.8 1];
    case 'contractility'
        ps = [0 0.05 0.1 0.15 0.3 0.45 0.6];
end
     
np = length(ps);
bins = 13;
gof_d = zeros(bins,np);
for pind=1:np
    fname = ['stats\stats_p-',param_name,'=', num2str(ps(pind)),'.mat'];
    if ~exist(fname, 'file')
        disp(["not found: ",fname]);
        continue;
    end
    load(fname,'gof');
    
    nc = length(gof);
    for i=0:bins-1
    	gof_d(i+1,pind) = sum(gof>=i & gof<i+1)*100/nc;
    end

end

heatmap(flipud(gof_d),'GridVisible','off','CellLabelColor','none');
ax = gca;
ax.XData = ps;
ax.YData = fliplr(1:bins);
colormap default
end