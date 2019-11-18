global timer

% lt: lower T2 threshold
% ha: higher A0 for HCs
% lr: lower HCs repulsion
% lrha: both the above
ext = '';
datestr = '131218';
wb = waitbar(0,'lower thresh');
for i=1:30
load(['stage1_complete-',datestr,'\lat',ext,'(',num2str(i) ,')_step(30).mat'],'g');
g = stage2differentiate(g);

g.area_feedback = 0;
g.fa = [1 1 1 1.2 1];
g = redistributeAreas(g);

timer = 0;
% parameters: area, tension, perimeter, mechanical feedback, HCs vdv, PCs
% pull, axial flow of SCs
g.paras = [4 ; 0.4; 0.1; 0; 1; 0; 0];
T1prob = 1;
T1eps = 0.1;
for t =1:101
    g = findTransitions(g,T1eps,T1prob,0.02); %0.05, 0.02, 0.02
    g = updateParameters(g);
    g = relaxLattice(g,5); %50

    g = killSmallCells(g,0.30);
    g = updateParameters(g);
    g = relaxLattice(g,5); %200
    
    disp(num2str(t))
    
    if mod(t-1,10) == 0
        save(['stage2_complete-',datestr,'/lat',ext,'(', num2str(i),')_step(', num2str((t-1)/10) ,')_lt.mat'],'g');
    end
end

waitbar(i/30);
end
close(wb)

wb = waitbar(0,'higher HC area');
for i=1:30
load(['stage1_complete-',datestr,'\lat',ext,'(',num2str(i) ,')_step(30).mat'],'g');
g = stage2differentiate(g);

g.area_feedback = 0;
g.fa = [1 1 1.2 1.2 1];
g = redistributeAreas(g);

timer = 0;
% parameters: area, tension, perimeter, mechanical feedback, HCs vdv, PCs
% pull, axial flow of SCs
g.paras = [4 ; 0.4; 0.1; 0; 1; 0; 0];
T1prob = 1;
T1eps = 0.1;
for t =1:101
    g = findTransitions(g,T1eps,T1prob,0.02); %0.05, 0.02, 0.02
    g = updateParameters(g);
    g = relaxLattice(g,5); %50

    g = killSmallCells(g,0.20);
    g = updateParameters(g);
    g = relaxLattice(g,5); %200
    
    disp(num2str(t))
    
    if mod(t-1,10) == 0
        save(['stage2_complete-',datestr,'/lat',ext,'(', num2str(i),')_step(', num2str((t-1)/10) ,')_ha.mat'],'g');
    end
end

waitbar(i/30);
end
close(wb)



wb = waitbar(0,'lower repulsion');
for i=1:30
load(['stage1_complete-',datestr,'\lat',ext,'(',num2str(i) ,')_step(30).mat'],'g');
g = stage2differentiate(g);

g.area_feedback = 0;
g.fa = [1 1 1 1.2 1];
g = redistributeAreas(g);

timer = 0;
% parameters: area, tension, perimeter, mechanical feedback, HCs vdv, PCs
% pull, axial flow of SCs
g.paras = [4 ; 0.4; 0.1; 0; 0.1; 0; 0];
T1prob = 1;
T1eps = 0.1;
for t =1:101
    g = findTransitions(g,T1eps,T1prob,0.02); %0.05, 0.02, 0.02
    g = updateParameters(g);
    g = relaxLattice(g,5); %50

    g = killSmallCells(g,0.20);
    g = updateParameters(g);
    g = relaxLattice(g,5); %200
    
    disp(num2str(t))
    
    if mod(t-1,10) == 0
        save(['stage2_complete-',datestr,'/lat',ext,'(', num2str(i),')_step(', num2str((t-1)/10) ,')_lr.mat'],'g');
    end
end

waitbar(i/30);
end
close(wb)


wb = waitbar(0,'lower repulsion and larger HCs');
for i=1:30
load(['stage1_complete-',datestr,'\lat',ext,'(',num2str(i) ,')_step(30).mat'],'g');
g = stage2differentiate(g);

g.area_feedback = 0;
g.fa = [1 1 1.2 1.2 1];
g = redistributeAreas(g);

timer = 0;
% parameters: area, tension, perimeter, mechanical feedback, HCs vdv, PCs
% pull, axial flow of SCs
g.paras = [4 ; 0.4; 0.1; 0; 0.1; 0; 0];
T1prob = 1;
T1eps = 0.1;
for t =1:101
    g = findTransitions(g,T1eps,T1prob,0.02); %0.05, 0.02, 0.02
    g = updateParameters(g);
    g = relaxLattice(g,5); %50

    g = killSmallCells(g,0.20);
    g = updateParameters(g);
    g = relaxLattice(g,5); %200
    
    disp(num2str(t))
    
    if mod(t-1,10) == 0
        save(['stage2_complete-',datestr,'/lat',ext,'(', num2str(i),')_step(', num2str((t-1)/10) ,')_lrha.mat'],'g');
    end
end

waitbar(i/30);
end
close(wb)

%% stats:

% plot mean of #of neighbors and chi^2 for hex fit
ts = 11; %timesteps
stg = 2; %stage 1 or 2

for si=1:4
    try
    switch si
        case 1
            add='lt';
        case 2
            add='ha';
        case 3
            add='lr';
        case 4
            add='lrha';
    end

non_t = zeros(ts,1);
gof_t = zeros(ts,1);
ar_t = zeros(ts,1); % area ratio of HC:SC

wb = waitbar(0, 'Winter is coming');
for t=0:ts-1
    non = [];
    gof = [];
    ar = [];
    for i=1:30
        load(['stage', num2str(stg),'_complete-',datestr,'\lat',ext,'(', num2str(i),')_step(', num2str(t),')_',add,'.mat'],'g');
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
            non = [non; ncount];
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
            chi2 = fit_to_hex(c_cent, ns_pos);
            gof_i = [gof_i; chi2];
        end
        gof = [gof; gof_i];
        
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
         ar = [ar; mean(ar_i)];
         
    end
    save(['stage', num2str(stg),'_complete-',datestr,'\stats',ext,'_step(', num2str(t),')_',add,'.mat'],'gof','non','ar');
    non_t(t+1) = mean(non);
    gof_t(t+1) = median(gof);
    ar_t(t+1) = mean(ar);
    waitbar((t+1)/ts);
end
close(wb);
save(['stage', num2str(stg),'_complete-',datestr,'\stats_all',ext,'_',add,'.mat'],'gof_t','non_t','ar_t');

t = 0:ts-1;
figure(10);
scatter(t, non_t);
xlabel('time (a.u.)');
ylabel('Number of SC neighbors');
set(gca, 'FontSize', 20, 'Box', 'on');
saveas(gcf,['plots\stage2_',add,'_non',ext],'png');
close

figure(20);
scatter(t, gof_t);
xlabel('time (a.u.)');
ylabel('Goodness of fit to hexagon');
set(gca, 'FontSize', 20, 'Box', 'on');
saveas(gcf,['plots\stage2_',add,'_hex',ext],'png');
close

figure(30);
scatter(t, ar_t);
xlabel('time (a.u.)');
ylabel('Ratio of HC:SC area');
set(gca, 'FontSize', 20, 'Box', 'on');
saveas(gcf,['plots\stage2_',add,'_area',ext],'png');
close
    catch
        disp('error');
    end
end
