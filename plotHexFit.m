
load(['D:\Users\TelAvivU-Analysis2\Desktop\Roie\vertex code\stage1_complete-131218\lat(10)_step(0).mat']);
g = stage2differentiate(g);
HCs = find(g.populations == 3);
% calc hex fit
gof = [];
bot = 0;
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
    [chi2, chi_err, fitted_params] = fit_to_hex(c_cent, ns_pos);
    gof = [gof; chi2, chi_err];
    
    continue; 
    
    % plot fit 
    LatticePresentation(g,0);
    hold on
    scatter(c_cent(1),c_cent(2));
    scatter(ns_pos(:,1), ns_pos(:,2));
    xlabel(['chi2: ', num2str(chi2)]);
    
    r = fitted_params(1);
    phi = fitted_params(2);
    alpha = fitted_params(3);
    
    % hexagon with the fitted parameters
    fitp = [ r*cos(0.5*phi),  r*sin(0.5*phi); ...
        r*cos(0.5*phi), -r*sin(0.5*phi); ...
        -r*cos(0.5*phi),  r*sin(0.5*phi); ...
        -r*cos(0.5*phi), -r*sin(0.5*phi); ...
        0             ,  2*r*sin(0.5*phi); ...
        0             , -2*r*sin(0.5*phi) ];
    
    % now rotate in 'alpha'
    R = [cos(alpha), sin(alpha); -sin(alpha) cos(alpha)];
    fitp = fitp*R;
    
    scatter(fitp(:,1) + c_cent(1), fitp(:,2) + c_cent(2), 'filled');  % fitted hexagon vertices 
    
    if bot ~= 27, [~, ~, bot] = ginput(1); end
    clf;
end
close