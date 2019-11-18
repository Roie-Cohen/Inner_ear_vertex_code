nlats = 20; % number of lattices
param_name = 'shear';
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

T2p = zeros(length(ps),1);
for pind=1:length(ps)
    T2tot = 0;
    for lat_i=1:nlats
        load(['stage1_complete-231218/lat(', num2str(lat_i),')_step(0).mat'],'g');
        gi=g;
        fname = ['stage1_', param_name,'\lat(', num2str(lat_i),')_p-',param_name,'=', num2str(ps(pind)),'_final.mat'];
        if ~exist(fname, 'file')
            disp(["not found: ",fname]);
            continue;
        end
        load(fname,'g');
        gf = g;
        T2s = sum(gf.dead - gi.dead);
        %     disp(['#T2 in lat ', num2str(lat_i), ' param ', num2str(ps(pind)), ' :', num2str(T2s)]);
        T2tot = T2tot + T2s;
    end
    T2p(pind) = T2tot;
    disp(['#T2 in simulations with ', param_name, ' = ', num2str(ps(pind)), ' :', num2str(T2tot)]);
end