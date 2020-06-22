%% create initial lattices
%% run this section to create random lattices
num_lat = 20;
for i=1:num_lat
    g = create_random_lattice(12,12);
    save(['random_lattices\rand_lat(',num2str(i),').mat']);
end

%% run this section to select initial pillar cells and inner hair cells