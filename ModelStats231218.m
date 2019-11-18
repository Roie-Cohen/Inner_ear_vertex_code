% plot mean of #of neighbors and chi^2 for hex fit
s1ts = 26; % numer of timesteps in stage 1
s2ts = 11; % numer of timesteps in stage 2
ts = s1ts + s2ts - 1; %timesteps
stg = 1; %stage 1 or 2
nlats = 50; % number of lattices
datestr = '271218';
ext = '';

non_t = zeros(ts,1);
gof_t = zeros(ts,1);
gof_t_med = zeros(ts,1);
ar_t = zeros(ts,1); % area ratio of HC:SC


wb = waitbar(0, 'Winter is coming');
for t=1:s1ts-1
    if t==s1ts
        stg = 2;
    end
    non = zeros(10000); non_ind = 1;
    gof = zeros(10000); gof_ind = 1;
    ar = zeros(10000);  ar_ind = 1;
    for i=1:nlats
        fname = ['stage', num2str(stg),'_complete-',datestr,'\lat(', num2str(i),')_step(', num2str(t-(s1ts-1)*(stg==2)),').mat'];
        if ~exist(fname, 'file'), continue; end
        load(fname,'g');

        % differentiate
        if stg==1, g = stage2differentiate(g); end
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
    ar = ar(ar~=0);
    save(['stats-',datestr,'\s`tats_step(', num2str(t),').mat'],'gof','non','ar');
    non_t(t+1) = mean(non);
    gof_t(t+1) = mean(gof);
    gof_t_med(t+1) = median(gof);
    ar_t(t+1) = mean(ar);
    waitbar((t+1)/ts);
end
close(wb);
% save(['stats-',datestr,'\stats_all.mat'],'gof_t','gof_t_med','non_t','ar_t');

%% plot
t = 0:ts-1;
figure(10);
scatter(t, non_t);
xlabel('time (a.u.)');
ylabel('Number of SC neighbors');
set(gca, 'FontSize', 20, 'Box', 'on');
figure(20);
scatter(t, gof_t_med);
xlabel('time (a.u.)');
ylabel('Goodness of fit to hexagon');
set(gca, 'FontSize', 20, 'Box', 'on');
figure(30);
scatter(t, ar_t);
xlabel('time (a.u.)');
ylabel('Ratio of HC:SC area');
set(gca, 'FontSize', 20, 'Box', 'on');
