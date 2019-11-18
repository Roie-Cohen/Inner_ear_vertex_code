function write_simulation_video_alternative(lat_num)

    time_steps = 100;
    sim_folder = 'alternative_2019-05-22';
    
    vid = VideoWriter(strcat('Simulations_vids/alternative_lat(',num2str(lat_num),').avi'));
    vid.Quality = 100; % 100 for best
    vid.FrameRate = 12;
    open(vid);
    figure(10);
    for t=0:time_steps
        load([sim_folder,'\lat(',num2str(lat_num),')_step(',num2str(t),').mat'],'g');
        LatticePresentation(g,0, 10);
        lm = 3.7;
        xlim([-lm lm]);
        ylim([-lm lm]);
        frame = getframe;
        writeVideo(vid,frame);
    end
    close(vid);
end