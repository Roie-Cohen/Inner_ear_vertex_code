
nlats = 20;
pind = 4;
param_name = 'tension';
ps = [0.05 0.1 0.2 0.3 0.45 0.6 0.8];

for lat_i=1:nlats
    fname = ['stage1_', param_name,'\lat(', num2str(lat_i),')_p-',param_name,'=', num2str(ps(pind)),'_final.mat'];
    if ~exist(fname, 'file')
        disp(["not found: ",fname]);
        continue;
    end
    load(fname,'g');
    g = semi_auto_s2diff(g);
    save(['stage1_', param_name,'\lat(', num2str(lat_i),')_p-',param_name,'=', num2str(ps(pind)),'_final_diff.mat'],'g');
end