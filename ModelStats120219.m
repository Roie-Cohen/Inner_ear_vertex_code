% plot mean of #of neighbors and chi^2 for hex fit
%ts = 26; % numer of timesteps in stage 1
nlats = 20; % number of lattices
global param_name
% param_name = 'shear';
switch param_name
    case 'area'
        mp = 1;
    case 'repulsion'
        mp = 5;
    case 'tension'
        mp = 2;
    case 'shear'
        mp = 7;
    case 'contractility'
        mp = 3;
end

my_ps = [6 ; 0.3; 0.15; 0; 1; 0.5; 0.4];
ps = (0.5:0.25:1.5)*my_ps(mp);

np = length(ps);

non_p = zeros(np,1);
gof_p = zeros(np,1);
gof_p_med = zeros(np,1);
ar_p = zeros(np,1); % area ratio of HC:SC

wb = waitbar(0, 'Winter is coming');
for pind=1:np 
    non = zeros(10000); non_ind = 1;
    gof = zeros(10000); gof_ind = 1;
    ar = zeros(10000);  ar_ind = 1;
    for i=1:nlats
        fname = ['stage1_', param_name,'\lat(', num2str(i),')_p-',param_name,'=', num2str(ps(pind)),'_final.mat'];
        if ~exist(fname, 'file')
            disp(["not found: ",fname]);
            continue;
        end
        load(fname,'g');

        % differentiate
        g = stage2differentiate(g);
        HCs = find(g.populations == 3);
        
        % calc #of neighbors
        for c=1:length(HCs)
            npop = g.populations(g.bonds(g.cells{HCs(c)+1}, 4));
            if ~isempty(find(npop~=2,1)), continue; end
            Lc = cellPerimeter(g, HCs(c));
            vidx = g.bonds(g.cells{HCs(c)+1},1);
            vert= getRelativePosition(g,vidx,HCs(c));
            nb=length(vidx);
            ncount = 0;
            for j=1:nb
                next = mod(j,nb)+1;   % the next vertex
                n2 = norm(vert(j,:)-vert(next,:));
                if (n2 > 0.03*Lc), ncount = ncount + 1; end
            end
            non(non_ind) = ncount;
            non_ind = non_ind + 1;
        end
        
        % calc hex fit
        gof_i = [];
        for c=1:length(HCs)
            npop = g.populations(g.bonds(g.cells{HCs(c)+1}, 4));
            if ~isempty(find(npop~=2,1)), continue; end
            sec_ns = [];
            for cn=1:length(HCs)
                if cn==c, continue; end
                cn_neigh = g.bonds(g.cells{HCs(cn)+1}, 4);
                if ~isempty(find(ismember(cn_neigh, g.bonds(g.cells{HCs(c)+1}, 4)),1))
                    sec_ns = [sec_ns; HCs(cn)];
                end
            end
            if length(sec_ns)< 4, continue; end
            ns_pos = zeros(length(sec_ns),2);
            c_cent = cellCenter(g, HCs(c));
            for j = 1:length(sec_ns)
                ns_pos(j, :) = cellCenter(g, sec_ns(j));
                if abs(ns_pos(j,1)-c_cent(1)) > 2*pi - abs(ns_pos(j,1)-c_cent(1)) % for periodic BC
                    ns_pos(j,1) = ns_pos(j,1) - 2*pi*sign(ns_pos(j,1));
                end
            end
            firstHCn = getFirstNeighbors(c_cent, ns_pos); % picking only first HC neighbors according to voronoi
            chi2 = fit_to_hex(c_cent, ns_pos(firstHCn, :));
            gof_i = [gof_i; chi2];
        end
        gof(gof_ind:gof_ind+length(gof_i)-1) = gof_i;
        gof_ind = gof_ind+length(gof_i);
        
        % calc HC:SC area ratio
        lat_areas = zeros(length(g.cells)-1);
        for c = 1:length(g.cells)-1
            lat_areas(c) = cellarea(g, c);
        end
        ar_i = [];
        for c=1:length(HCs)
            c_area = lat_areas(HCs(c));
            neigh = g.bonds(g.cells{HCs(c)+1}, 4);
            for ni = 1:length(neigh)
                if g.populations(neigh(ni))==2
                    ar_i = [ar_i; c_area/lat_areas(neigh(ni))];
                end
            end
        end
         ar(ar_ind) = mean(ar_i);
         ar_ind = ar_ind + 1;
         
    end
    non = non(non~=0);
    gof = gof(gof~=0);
    gof = gof(~isinf(gof));
    ar = ar(ar~=0);
    save(['stats\stats_p-',param_name,'=', num2str(ps(pind)),'.mat'],'gof','non','ar');
    non_p(pind) = mean(non);
    gof_p(pind) = mean(gof);
    gof_p_med(pind) = median(gof);
    ar_p(pind) = mean(ar);
    waitbar(pind/np);
end
close(wb);
save(['stats\stats_all_',param_name,'.mat'],'gof_p','gof_p_med','non_p','ar_p');

%% plot
switch param_name
    case 'area'
        x_name = '\alpha';
    case 'repulsion'
        x_name = '\beta';
    case 'tension'
        x_name = '\gamma';
    case 'shear'
        x_name = '\zeta';
    case 'contractility'
        x_name = '\Gamma';
end

figure(10);
scatter(ps, non_p);
xlabel(x_name);
ylabel('Number of SC neighbors');
set(gca, 'FontSize', 20, 'Box', 'on');
saveas(gcf,['plots\non(',param_name,').png']);
saveas(gcf,['plots\non(',param_name,').fig']);
% xlim([0,ts-1]);
% xticks(0:5:35)
figure(20);
scatter(ps, gof_p);
xlabel(x_name);
ylabel('Goodness of fit to hexagon');
set(gca, 'FontSize', 20, 'Box', 'on');
saveas(gcf,['plots\hex(',param_name,').png']);
saveas(gcf,['plots\hex(',param_name,').fig']);
% xlim([0,ts-1]);
% xticks(0:5:35)
figure(30);
scatter(ps, ar_p);
xlabel(x_name);
ylabel('Ratio of HC:SC area');
set(gca, 'FontSize', 20, 'Box', 'on');
saveas(gcf,['plots\ar(',param_name,').png']);
saveas(gcf,['plots\ar(',param_name,').fig']);
% xlim([0,ts-1]);
% xticks(0:5:35)