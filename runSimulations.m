%% run the simulation for the lattices in initial_lattices folder

num_lat = 20; % number of initial lattices to run the simulation on

wb = waitbar(0, 'Simulations run');
for i=1:num_lat
    load(['initial_lattices\init_lat(',num2str(i),').mat']);
    g = stage1(g, 400, i);
    g = stage2(g, 100, i);
    waitbar(i/num_lat);
end
close(wb);