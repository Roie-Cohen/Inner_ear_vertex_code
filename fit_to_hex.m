function [chi2, chi_err, fitted_params] = fit_to_hex(center, points)

% center the system to (0,0) 
 for i=1:2
     points(:,i) = points(:,i) - center(i);
 end

% hexagon parameters:
% r - distance from center to nearest vertex
% phi - opening angle of the hexagon (60 deg in regular hex)
% alpha - orientation of the main axis of the hexagon
np = size(points, 1);
cmin = 1E50;
dists = sqrt(sum(points.^2,2));
dmax = max(dists);
dmin = min(dists);
davg = mean(dists);
dstep = (dmax-dmin)/20;
for r=dmin:dstep:dmax
    for phi=pi/18:pi/36:17*pi/18
        ares = 36;
        astep = pi/ares;
        R = [cos(astep), sin(astep); -sin(astep) cos(astep)];
        fitp = [ r*cos(0.5*phi),  r*sin(0.5*phi); ...
            r*cos(0.5*phi), -r*sin(0.5*phi); ...
            -r*cos(0.5*phi),  r*sin(0.5*phi); ...
            -r*cos(0.5*phi), -r*sin(0.5*phi); ...
            0             ,  2*r*sin(0.5*phi); ...
            0             , -2*r*sin(0.5*phi) ];
        for a=1:ares
            fitp = fitp*R;
            sumd = 0;
            for j=1:np
                d2 = min(sum((fitp-points(j,:)).^2,2));
                sumd = sumd + d2;
            end
            sumd = sumd/((0.1*davg)^2);
            if sumd < cmin
                cmin = sumd;
                fitted_params = [r, phi, astep*a];
            end
        end
    end
end
% returning chi2 reduced
chi2 = cmin/(np-3);
chi_err = sqrt(2/(np-3));

% disp(['Var = ',num2str((0.1*davg)^2)]); %remove

 