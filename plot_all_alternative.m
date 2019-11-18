nlats = 20;
npar = 5;
nvar = 5;
gofstats = zeros(npar, nvar, 2);
nonstats = zeros(npar, nvar, 2);
for pari = 1:npar

switch pari
    case 1 
        param_name = 'area';
        mp = 1;
    case 2 
        param_name = 'repulsion';
        mp = 5;
    case 3 
        param_name = 'tension';
        mp = 2;
    case 4
        param_name = 'shear';
        mp = 7;
    case 5
        param_name = 'contractility';
        mp = 3;
end

my_ps = [6 ; 0.3; 0.15; 0; 1; 0.5; 0.4];
ps = (0.5:0.25:1.5)*my_ps(mp);


% plot gof
np = length(ps);

for pind=1:np
    fname = ['stats\stats_p-',param_name,'=', num2str(ps(pind)),'.mat'];
    if ~exist(fname, 'file')
        disp(["not found: ",fname]);
        continue;
    end
    load(fname,'gof');
    gofstats(pari, pind, :) = [mean(gof), std(gof)/sqrt(length(gof))];
%     color = 'b';
%     scatter(ps(pind), mean(gof),color,'filled');
%     errorbar(ps(pind), mean(gof), std(gof)/sqrt(length(gof)), color);
%     
end
% gofstats(pari, :, :) = gofstats(pari, :, :) / gofstats(pari, 3, 1); % normalize according to the original simulated value

% ylim([4 12]);
% xlim([ps(1) ps(length(ps))]);
% xlabel(x_name);
% ylabel('Goodness of fit to hexagon');
% set(gca, 'FontSize', 20, 'Box', 'on');
% saveas(gcf,['plots\hex(',param_name,')_scatter.png']);
% saveas(gcf,['plots\hex(',param_name,')_scatter.fig']);


% plot non
np = length(ps);
for pind=1:np
    fname = ['stats\stats_p-',param_name,'=', num2str(ps(pind)),'.mat'];
    if ~exist(fname, 'file')
        disp(["not found: ",fname]);
        continue;
    end
    load(fname,'non');
    nonstats(pari, pind, :) = [mean(non), std(non)/sqrt(length(non))];  
end
% nonstats(pari, :, :) = nonstats(pari, :, :) / nonstats(pari, 3, 1); % normalize according to the original simulated value

end

figure(10)
hold on
figure(20)
hold on
for pari = 1:npar
    switch pari
        case 1
            color = [0 0.45 0.75];
            dispName = 'Area';
        case 2
            color = [0.85 0.33 0.10];
            dispName = 'HCs repulsion';
        case 3
            color = [0.93 0.69 0.13];
            dispName = 'Junctional tension';
        case 4
            color = [0.49 0.18 0.56];
            dispName = 'Shear';
        case 5
            color = [0.47 0.67 0.19];
            dispName = 'Contractility';
    end
    
    figure(10);
%     scatter(0.5:0.25:1.5, gofstats(pari,:,1), 'MarkerFaceColor', color);
    errorbar(0.5:0.25:1.5, gofstats(pari,:,1), gofstats(pari,:,2),'-o','MarkerSize',8,...
        'MarkerFaceColor',color,'LineStyle','none','LineWidth',2,'CapSize', 15,'DisplayName',dispName);
    figure(20)
%     scatter(0.5:0.25:1.5, nonstats(pari,:,1), 'MarkerFaceColor', color);
    errorbar(0.5:0.25:1.5, nonstats(pari,:,1), nonstats(pari,:,2),'-o','MarkerSize',8,...
        'MarkerFaceColor',color,'LineStyle','none','LineWidth',2,'CapSize', 15,'DisplayName',dispName);
end
figure(10);
xlim([0.45 1.55]);
xticks([0.5:0.25:1.5]);
ylim([5 12]);
xlabel('p/p_0');
ylabel('Goodness of fit to hexagon');
set(gca, 'FontSize', 20, 'Box', 'on');
legend('FontSize',12)
figure(20);
xlim([0.45 1.55]);
xticks([0.5:0.25:1.5]);
ylim([4 6]);
xlabel('p/p_0');
ylabel('Number of SC neighbors');
set(gca, 'FontSize', 20, 'Box', 'on');
legend('FontSize',12)