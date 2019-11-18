function produce_stats(sim_folder, nsim, nsteps, stage)
% produces the statistics files for the simulations in the folder 'sim_folder'
% nsim - number of simulations. nsteps - number of time steps per simulation.

stat_every = 2; %10 take statistic every # steps
stat_steps = floor(nsteps/stat_every); % number of statistics points
for t=0:stat_steps
    
    gs = cell(nsim, 1); % will include all the lattices of timepoint t*stat_every
    for i=1:nsim
        fname = [sim_folder,'\lat(', num2str(i),')_step(', num2str(t*stat_every),').mat'];
        if ~exist(fname, 'file')
            disp(["not found: ",fname]);
            continue;
        end
        load(fname,'g');
        % addition of statistics fields to 'g'
        g.non = zeros(500); % number of SC neighbors for HCs
        g.gof = zeros(500); % goodness of fit to hexagon
        g.ar  = zeros(500); % area ratio of HCs to SCs
        g.arr  = zeros(500); % real area ratio of HCs to SCs (including curvature)
        gs{i} = g;
    end
    
    parfor i=1:nsim
        
        non_ind = 1;
        gof_ind = 1;
        ar_ind = 1;
        
        g = gs{i};
        % differentiate
        if stage == 1 % for stage1 this differentiation if essential to define the cells in the middel row
            g = stage2differentiate(g); 
        end
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
        if ~isnan(g.verts(1,1))
            real_areas= getRealAreas(g);
            lat_areas = zeros(length(g.cells)-1);
            for c = 1:length(g.cells)-1
                lat_areas(c) = cellarea(g, c);
            end
            g.centroid = cellCOM(g,1:length(g.cells)-1);
            ar_i = [];
            arr_i = [];
            for c=1:length(HCs)
                npop = g.populations(g.bonds(g.cells{HCs(c)+1}, 4));
                if ~isempty(find(npop~=2,1)), continue; end
                c_area = lat_areas(HCs(c));
                cr_area = real_areas(HCs(c));
                neigh = g.bonds(g.cells{HCs(c)+1}, 4);
                for ni = 1:length(neigh)
                    if g.populations(neigh(ni))==2
                        ar_i = [ar_i; c_area/lat_areas(neigh(ni))];
                        arr_i = [arr_i; cr_area/real_areas(neigh(ni))];
                    end
                end
                g.ar(ar_ind) = mean(ar_i);
            g.arr(ar_ind) = mean(arr_i);
            ar_ind = ar_ind + 1;
            end
        else
%              disp(num2str(i)); 
        end
        
        g.non = g.non(g.non~=0);
        g.gof = g.gof(g.gof~=0);
        g.gof = g.gof(~isinf(g.gof));
        g.ar = g.ar(g.ar~=0);
        g.ar = g.ar(~isnan(g.ar));
        g.arr = g.arr(g.arr~=0);
        g.arr = g.arr(~isnan(g.arr));
        gs{i} = g;
    end
    
    % make statistics folder if not exist
    stat_folder = [sim_folder, '_stats'];
    if ~exist(stat_folder, 'dir')
        mkdir(stat_folder);
    end
    % save statistics and make average statistics
    non = [];
    gof = [];
    ar = [];
    arr = [];
    for i=1:nsim
        g = gs{i};
        fname = [stat_folder,'\stat_lat(', num2str(i),')_step(', num2str(t),').mat'];
        save(fname,'g');
        
        non = [non; g.non];
        gof = [gof; g.gof];
        ar = [ar; g.ar];
        arr = [arr; g.arr];
    end
    
    non_t = [mean(non), std(non)/sqrt(length(non))];
    gof_t = [mean(gof), std(gof)/sqrt(length(gof))];
    ar_t = [mean(ar), std(ar)/sqrt(length(ar))];
    arr_t = [mean(arr), std(arr)/sqrt(length(arr))];
    
    % save average statistics for this timestep
    fname = [stat_folder,'\stat_step(', num2str(t),').mat'];
    save(fname,'non_t','gof_t','ar_t','arr_t');
    
    disp(['Producing statistics - step ', num2str(t), '. Percentage: ', num2str(round(100*(t+1)/(stat_steps+1) ,2))]);
end
    
end
