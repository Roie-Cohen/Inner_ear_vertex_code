%% make simulation movies

lat_inds = 1:20; % the lattices numbers to make a movie for
stg1_ts = 40;
stg2_ts = 10;

figure(100)

for i=lat_inds
    v = VideoWriter(['simulation_movies\lat(', num2str(i),').avi']);
    v.FrameRate = 5;
    
    open(v);
    for t=0:stg1_ts-1
        load(['stage1_timepoints\lat(', num2str(i), ')_step(',num2str(t),').mat']);
        clf;
        LatticePresentation(g, 0);
        frm = getframe(gcf);
        writeVideo(v,frm);
    end
    for t=1:stg2_ts
        load(['stage2_timepoints\lat(', num2str(i), ')_step(',num2str(t),').mat']);
        clf;
        LatticePresentation(g, 0);
        frm = getframe(gcf);
        writeVideo(v,frm);
    end
    close(v);
    clear v;
end

close
