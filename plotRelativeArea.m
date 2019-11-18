% ablated_cells = [4 125 5 104 126 115 139 131 142 119 121 109 2 133];
stage2_fold = 'stage2_2019-05-22'; % folder with stage1 simulations
s2tf = 100; % final step in stage1
HCs_a = []; SCs_a = [];
for lat_i = 1:50
    
    lat_file_name = [stage2_fold,'/lat(', num2str(lat_i),')_step(',num2str(s2tf),').mat'];
    if exist(lat_file_name, 'file')
        load([stage2_fold,'/lat(', num2str(lat_i),')_step(',num2str(s2tf),').mat'],'g'); % initial
    else
        continue;
    end
    gi = g; % initial lattice
    
    
    for c = 1:144
        ab_file_name = ['ablated_cells/ablation-lat(',num2str(lat_i),')-cell(',num2str(c),').mat'];
        if exist(ab_file_name, 'file')
            load(ab_file_name,'g');
            gf = g;
        else
            continue;
        end
        
        
        % SC and HC first neighbors of ablated cell
        n_bonds = gi.cells{c+1};
        ns = gi.bonds(n_bonds, 4); % first neighbors
        %     second_ns = first_ns;
        %     for j=1:length(first_ns)
        %         cc = first_ns(j);
        %         second_ns = [second_ns; g.bonds(g.cells{cc+1}, 4)];
        %     end
        %     ns = unique(second_ns);
        types = gi.populations(ns);
        ns = ns(types == 3 | types == 2);
        ns = ns(ns~=c);
        types = gi.populations(ns); % matching types (SC or HC) in 'ns'
        
        % initial area of the surrounding cells
        areas_i = zeros(length(ns), 1);
        for j=1:length(ns)
            areas_i(j) = cellarea(gi, ns(j));
        end
        
        % final area
        areas_f = zeros(length(ns), 1);
        for j=1:length(ns)
            areas_f(j) = cellarea(gf, ns(j));
        end
        
        % relative areas of HCs and SCs
        ar = areas_f./areas_i;
        a_diff = abs(1 - ar);
        HCs_a = [HCs_a; a_diff(types == 3)];
        SCs_a = [SCs_a; a_diff(types == 2)];
        
    end
end
HCs_a = HCs_a(~isnan(HCs_a));
SCs_a = SCs_a(~isnan(SCs_a));

figure(100);
hold on
ylabel('|1-A/A0|');
set(gca, 'FontSize', 24);
set(gca, 'Box', 'on');
X = categorical({'HCs','SCs'});
X = reordercats(X,{'HCs','SCs'});
Y = [mean(HCs_a), mean(SCs_a)];
bar(X,Y)
Yerr = [std(HCs_a)/sqrt(length(HCs_a)) std(SCs_a)/sqrt(length(SCs_a))];
er = errorbar(Y,Yerr);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 2;

[~,pval,~,~] = ttest2(HCs_a,SCs_a);
disp(['two-sample t-test p-value: ', num2str(pval)]);

[~,pval] = kstest2(HCs_a,SCs_a);
disp(['two-sample ks-test p-value: ', num2str(pval)]);