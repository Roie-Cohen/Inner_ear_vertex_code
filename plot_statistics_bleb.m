ts = 11; %timesteps
datestr = '2019-05-22';

non = zeros(ts,2);
gof = zeros(ts,2);
ar = zeros(ts,2); % area ratio of HC:SC

for t=0:ts-1
   
    load(['bleb_',datestr,'_stats\stat_step(', num2str(t),').mat'],'gof_t','non_t','ar_t');
    non(t+1,:) = non_t;
    gof(t+1,:) = gof_t;
    ar(t+1,:) = ar_t;
    
end

% save(['stats-',datestr,'\stats_all.mat'],'gof_t','gof_t_med','non_t','ar_t');

%% plot 
t = 0:ts-1;
figure(10);
errorbar(t,non(:,1),non(:,2),'o','linewidth',2);
xlabel('time (a.u.)');
ylabel('Number of SC neighbors');
set(gca, 'FontSize', 20, 'Box', 'on');
ylim([4.5 6]);
xlim([0 ts-1]);

figure(20);
errorbar(t,gof(:,1),gof(:,2),'o','linewidth',2);
xlabel('time (a.u.)');
ylabel('Goodness of fit to hexagon');
set(gca, 'FontSize', 20, 'Box', 'on');
ylim([3 12]);
xlim([0 ts-1]);

figure(30);
errorbar(t,ar(:,1),ar(:,2),'o','linewidth',2);
xlabel('time (a.u.)');
ylabel('Ratio of HC:SC area');
set(gca, 'FontSize', 20, 'Box', 'on');
ylim([0.5 2])