function kappa = getCurvature(g, bo)
% returns the curvature of the bond 'bo' (calculated with laplace law)

alpha_g = g.paras(1); % global factor for compressibility
gamma_g = g.paras(2); % global factor for tension

gamma_nm = g.bonds(bo, 5); % tension of the bond

cs = g.bonds(bo, 3:4);     % neighboring cells' indices
n = cs(1); m = cs(2);      % cells' indices
alpha_n = g.cmpr(n);   % compressibility of cell n
alpha_m = g.cmpr(m);   % compressibility of cell m
dAn = cellarea(g, n) - g.areas(n); % An-A0n
dAm = cellarea(g, m) - g.areas(m); % Am-A0m

kappa = (alpha_g/gamma_g)*(alpha_n*dAn - alpha_m*dAm)/gamma_nm; % curvature

end