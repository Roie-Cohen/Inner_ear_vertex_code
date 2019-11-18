ri = 0; rf = 4;
r = linspace(ri,rf,1000);
close all
figure;
hold on
for l=1:6
    k=2*l;
    V = @(r) (sigma./r).^k - (sigma./r).^l;
    plot(r,V(r))
end
ylim([-0.5,0.5])
