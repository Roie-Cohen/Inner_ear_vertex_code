function write_simulation_video_params(lat_num, timesteps, param_name, param_val)

    sim_folder = ['stage1_2019-05-22_',param_name];
    
    vid = VideoWriter(strcat('Simulations_vids/lat(',num2str(lat_num),'_',param_name,'=',num2str(param_val),'.avi'));
    vid.Quality = 50; % 100 for best
    open(vid);
    figure(11);
    for t=0:timesteps
        load([sim_folder,'\lat(',num2str(lat_num),')_p=',num2str(param_val),'_step(',num2str(t),').mat'],'g');
        LatticePresentation(g,0, 11);
        lm = 3.7;
        xlim([-lm lm]);
        ylim([-lm lm]);
        frame = getframe;
        writeVideo(vid,frame);
    end
    close(vid);
end