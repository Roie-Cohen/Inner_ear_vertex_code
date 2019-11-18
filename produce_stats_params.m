function produce_stats_params(sim_folder, nsim, last_step, param_name, param_values, save_folder)
% produces the statistics files for the simulations with different
% parameters values for the last step in stage1.
% nsim - number of simulations. last_step - the last step of the
% simulations. param_name - parameter name as indicated in runStage1.m
% param_values - list of the values 

for pari = 1:length(param_values)
    gs = cell(nsim, 1); % will include all the lattices in the last timepoint of the simulation
    p = param_values(pari);
    for i=1:nsim
        fname = [sim_folder,'\lat(', num2str(i),')_p=', num2str(p),'_step(', num2str(last_step),').mat'];
        if ~exist(fname, 'file')
            disp(["not found: ",fname]);
            continue;
        end
        load(fname,'g');
        % addition of statistics fields to 'g'
        g.non = zeros(500); % number of SC neighbors for HCs
        g.gof = zeros(500); % goodness of fit to hexagon
        g.ar  = zeros(500); % area ratio of HCs to SCs
        g.damaged = 0; % indicates if there's a problem with the lattice
        gs{i} = g;
    end
    
    parfor i=1:nsim
        
        non_ind = 1;
        gof_ind = 1;
        ar_ind = 1;
        
        g = gs{i};
        
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
            g.non(non_ind) = ncount;
            non_ind = non_ind + 1;
        end
        
        % calc hex fit
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
            c_cent = cellCOM(g, HCs(c));
            if isnan(c_cent)
                disp(['problem with simulation number ', num2str(i), ' and parameter value ', num2str(p)]);
                g.damaged = 1;
                continue;
            end
            for j = 1:length(sec_ns)
                ns_pos(j, :) = cellCOM(g, sec_ns(j));
                if abs(ns_pos(j,1)-c_cent(1)) > 2*pi - abs(ns_pos(j,1)-c_cent(1)) % for periodic BC
                    ns_pos(j,1) = ns_pos(j,1) - 2*pi*sign(ns_pos(j,1));
                end
            end
            firstHCn = getFirstNeighbors(c_cent, ns_pos); % picking only first HC neighbors according to voronoi
            chi2 = fit_to_hex(c_cent, ns_pos(firstHCn, :));
            g.gof(gof_ind) = chi2;
            gof_ind = gof_ind+1;
        end

        
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
            g.ar(ar_ind) = mean(ar_i);
            ar_ind = ar_ind + 1;
        end
        
        g.non = g.non(g.non~=0);
        g.gof = g.gof(g.gof~=0);
        g.gof = g.gof(~isinf(g.gof));
        g.ar = g.ar(g.ar~=0);
        g.ar = g.ar(~isnan(g.ar));
        gs{i} = g;
    end
    
    % make statistics folder if not exist
    if ~exist(save_folder, 'dir')
        mkdir(save_folder);
    end
    % save statistics and make average statistics
    non = [];
    gof = [];
    ar = [];
    for i=1:nsim
        g = gs{i};
        fname = [save_folder,'\stat_lat(', num2str(i),')_',param_name,'=',num2str(p),'_final.mat'];
        save(fname,'g');
        
        if ~g.damaged
            non = [non; g.non];
            gof = [gof; g.gof];
            ar = [ar; g.ar];
        end
    end
    
    non_t = [mean(non), std(non)/sqrt(length(non))];
    gof_t = [mean(gof), std(gof)/sqrt(length(gof))];
    ar_t = [mean(ar), std(ar)/sqrt(length(ar))];
    
    % save average statistics
    fname = [save_folder,'\stat_',param_name,'=',num2str(p),'.mat'];
    save(fname,'non_t','gof_t','ar_t');
    
    disp(['Producing statistics - Percentage: ', num2str(100*pari/length(param_values))]);
   
end
end
