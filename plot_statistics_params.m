stat_folder = 'stage1_2019-05-22_params-stats';

param_names = ["gl-alpha","gl-gamma","gl-Gama","gl-eta","gl-zeta","gl-sigma","HC-alpha","HC-Gamma"];
param_signs = ["\alpha_0","\gamma_0","\Gamma_0","\eta","\zeta","\sigma","\alpha_H_C","\Gamma_H_C"];
p0 = [12, 0.3, 0.15, 0.5, 1, 1, 2, 2];
param_vals = 0.5:0.25:1.5;
pn = length(param_vals);

for i=1:length(param_names)
non = zeros(pn,2); 
gof = zeros(pn,2);
ar = zeros(pn,2); % area ratio of HC:SC

for pind=1:pn
    load([stat_folder, '\stat_', char(param_names(i)),'=', num2str(param_vals(pind)*p0(i)),'.mat'],'gof_t','non_t','ar_t');
    non(pind,:) = non_t;
    gof(pind,:) = gof_t;
    ar(pind,:) = ar_t;
    
end
 
figure(10);
hold on
errorbar(param_vals,non(:,1),non(:,2),'o','linewidth',2,'DisplayName',param_signs(i));
xlabel('p/p_0');
ylabel('Number of SC neighbors');
set(gca, 'FontSize', 20, 'Box', 'on');
ylim([4 6]);

figure(20);
hold on
errorbar(param_vals,gof(:,1),gof(:,2),'o','linewidth',2,'DisplayName',param_signs(i));
xlabel('p/p_0');
ylabel('Goodness of fit to hexagon');
set(gca, 'FontSize', 20, 'Box', 'on');
ylim([3 12]);

figure(30);
hold on
errorbar(param_vals,ar(:,1),ar(:,2),'o','linewidth',2,'DisplayName',param_signs(i));
xlabel('p/p_0');
ylabel('Ratio of HC:SC area');
set(gca, 'FontSize', 20, 'Box', 'on');
ylim([0 2])

end

% for i=10:10:30
%     figure(i);
%     legend(param_signs(1),param_signs(2),param_signs(3),param_signs(4),param_signs(5),param_signs(6),param_signs(7),param_signs(8));
% end