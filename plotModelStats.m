

%% 
% plotting number of SCs neighbors of middle row HCs.
figure(10)
non = [];
for i=1:30
    load(['stage1_complete-211118\lat(', num2str(i),')_step(20).mat'],'g');
    g = stage2differentiate(g);
    HCs = find(g.populations == 3);
    for c=1:length(HCs)
        npop = g.populations(g.bonds(g.cells{HCs(c)+1}, 4));
        if ~isempty(find(npop~=2,1)), continue; end
        non = [non; length(npop)];
    end
    disp(num2str(i));
end
hist(non);
sim_non = non;
hold on
%% 
% plotting number of neighbors for arbitrary distribution of lateral
% inhibition.

non = [];
load('lattices/LIlat_autodiff1','g');
gorig = g;
for i=1:50
    g = gorig;
    g = LIdifferentiation(g, 2.4);
    HCs = find(g.populations == 3);
    for c=1:length(HCs)
        npop = g.populations(g.bonds(g.cells{HCs(c)+1}, 4));
        if ~isempty(find(npop~=2,1)), continue; end
        non = [non; length(npop)];
    end
    disp(num2str(i));
end
hist(non);

%% 
% plotting the goodness of fit to hexagon.

figure(20)
gof = [];
for i=1:30
    gof_i = [];
    load(['stage1_complete-211118\lat(', num2str(i),')_step(20).mat'],'g');
    g = stage2differentiate(g);
    HCs = find(g.populations == 3);
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
        chi2 = fit_to_hex(c_cent, ns_pos);
        gof_i = [gof_i; chi2];
    end
    gof = [gof; gof_i];
    disp([num2str(i),': ',num2str(mean(gof_i))]);
end
sim_gof = gof;
hist(gof);
hold on

%% 
% plotting the goodness of fit to hexagon for arbitrary distribution of lateral
% inhibition.

gof = [];
load('lattices/LIlat_autodiff1','g');
gorig = g;
for i=1:50
    gof_i = [];
    g = gorig;
    g = LIdifferentiation(g, 2.4);
    g = stage2differentiate(g);
    HCs = find(g.populations == 3);
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
        chi2 = fit_to_hex(c_cent, ns_pos);
        gof_i = [gof_i; chi2];
    end
    gof = [gof; gof_i];
    disp([num2str(i),': ',num2str(mean(gof_i))]);
end
hist(gof);