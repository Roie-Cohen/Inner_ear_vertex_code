s1ts = 41; %26 numer of timesteps in stage 1
s2ts = 11; % numer of timesteps in stage 2. s2ts=1 for stage1 only
ts = s1ts + s2ts - 1; %timesteps
datestr = '2019-05-22';

non = zeros(ts,2);
gof = zeros(ts,2);
ar = zeros(ts,2); % area ratio of HC:SC
arr = zeros(ts,2); % real area ratio of HC:SC (including curvature)

stg = 1; %stage 1 or 2
for t=0:ts-1
    if t==s1ts
        stg = 2;
    end
   
    load(['stage', num2str(stg),'_',datestr,'_stats\stat_step(', num2str(t-(s1ts-1)*(stg==2)),').mat'],'gof_t','non_t','ar_t','arr_t');
    non(t+1,:) = non_t;
    gof(t+1,:) = gof_t;
    ar(t+1,:) = ar_t;
    arr(t+1,:) = arr_t;
    
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

% figure(30);
% % size at the end of the stages
% c = categorical({'Stage1','Stage2','Stage2*'});
% arf = [ar(s1ts,1), ar(ts,1), arr(ts,1)];
% arf_err = [ar(s1ts,2), ar(ts,2), arr(ts,2)];
% bar(c,arf)
% hold on
% errorbar(arf,arf_err,'o','linewidth',2, 'color', [0.15, 0.15, 0.15], 'CapSize', 12);
% xlabel('time (a.u.)');
% ylabel('Ratio of HC:SC area');
% set(gca, 'FontSize', 20, 'Box', 'on');
% ylim([0 4])

figure(30);
% size at the end of the stages
c = categorical({'Stage1','Stage2'});
arf = [ar(s1ts,1), ar(ts,1)];
arf_err = [ar(s1ts,2), ar(ts,2)];
bar(c,arf)
hold on
errorbar(arf,arf_err,'o','linewidth',2, 'color', [0.15, 0.15, 0.15], 'CapSize', 12);
xlabel('time (a.u.)');
ylabel('Ratio of HC:SC area');
set(gca, 'FontSize', 20, 'Box', 'on');
ylim([0.5 1.3])

figure(40);
% size at the end of the stages
c = categorical({'Stage1','Stage2'});
arrf = [arr(s1ts,1), arr(ts,1)];
arrf_err = [arr(s1ts,2), arr(ts,2)];
bar(c,arrf)
hold on
errorbar(arrf,arrf_err,'o','linewidth',2, 'color', [0.15, 0.15, 0.15], 'CapSize', 12);
xlabel('time (a.u.)');
ylabel('Ratio of HC:SC area');
set(gca, 'FontSize', 20, 'Box', 'on');
ylim([0 4])