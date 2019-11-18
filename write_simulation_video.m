function write_simulation_video(sim_num, time_steps, stg)

    sim_folder = ['stage', num2str(stg),'_2019-05-22'];
    
    vid = VideoWriter(strcat('Simulations_vids/stage',num2str(stg),'-lat(',num2str(sim_num),').avi'));
    vid.Quality = 100; % 100 for best
    vid.FrameRate = 12;
    open(vid);
    figure(10);
    for t=0:time_steps
        load([sim_folder,'\lat(',num2str(sim_num),')_step(',num2str(t),').mat'],'g');
        LatticePresentation(g,0, 10);
        lm = 3.7;
        xlim([-lm lm]);
        ylim([-lm lm]);
        frame = getframe;
        writeVideo(vid,frame);
    end
    close(vid);
end