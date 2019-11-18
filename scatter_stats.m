function scatter_stats(param_name)

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

nlats = 20;
T2thresh = 5;
% count T2
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
        T2tot = T2tot + T2s;
    end
    T2p(pind) = T2tot;
end


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

% plot gof
ds = min(diff(ps))*0.5;
np = length(ps);
figure(10);
hold on
for pind=1:np
    fname = ['stats\stats_p-',param_name,'=', num2str(ps(pind)),'.mat'];
    if ~exist(fname, 'file')
        disp(["not found: ",fname]);
        continue;
    end
    load(fname,'gof');
%     scatter(ps(pind)+ds*(rand(length(gof),1)-0.5), gof,'b+');
%     scatter(ps(pind), mean(gof),'b');
%     if T2p(pind) >= T2thresh
%         color = 'r';
%     else
%         color = 'k';
%     end
    color = 'b';
    scatter(ps(pind), mean(gof),color,'filled');
    errorbar(ps(pind), mean(gof), std(gof)/sqrt(length(gof)), color);
    
end
ylim([4 12]);
xlim([ps(1) ps(length(ps))]);
xlabel(x_name);
ylabel('Goodness of fit to hexagon');
set(gca, 'FontSize', 20, 'Box', 'on');
saveas(gcf,['plots\hex(',param_name,')_scatter.png']);
saveas(gcf,['plots\hex(',param_name,')_scatter.fig']);


% plot non
np = length(ps);
figure(20);
hold on
for pind=1:np
    fname = ['stats\stats_p-',param_name,'=', num2str(ps(pind)),'.mat'];
    if ~exist(fname, 'file')
        disp(["not found: ",fname]);
        continue;
    end
    load(fname,'non');
%     scatter(ps(pind)+ds*(rand(length(gof),1)-0.5), gof,'b+');
%     scatter(ps(pind), median(non),'k');
%     if T2p(pind) >= T2thresh
%         color = 'r';
%     else
%         color = 'b';
%     end
    color = 'b';
    scatter(ps(pind), mean(non),color,'filled');
    errorbar(ps(pind), mean(non), std(non)/sqrt(length(non)), color);
    
end
ylim([4 6]);
xlim([ps(1) ps(length(ps))]);
xlabel(x_name);
ylabel('Number of SC neighbors');
set(gca, 'FontSize', 20, 'Box', 'on');
saveas(gcf,['plots\non(',param_name,')_scatter.png']);
saveas(gcf,['plots\non(',param_name,')_scatter.fig']);
end