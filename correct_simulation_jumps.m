%% initial parameters
% correct jumps in simulation movies
% frames g structure files are in this folder:
fold = 'E:\Roie\vertex_code\stage1_2019-05-22';
% the format of file names are
name_for = 'lat(lll)_step(sss).mat';
% where 'lll' replaces by latice/simulation number and 'sss' replaced by
% frame/step number.

sim_num = 23;
name_for = strrep(name_for, 'lll', num2str(sim_num));
sim_length =400;

% the code follows the movement of cell 'c' and marks the frames where
% there was a large jump in the positions of 'c'.
c = 127; % we can change for any other cell
jthresh = 0.2; % jump threshold. 
marked_frames = []; 
gs = cell(sim_length + 1,1);

%% find jumps
i = 0;
filename = [fold,'\',strrep(name_for, 'sss', num2str(i))];
load(filename,'g');
for k=1:length(g.cells)-1
    if ~g.dead(k), g.centroid(k,:) = cellCOM(g,k); end
end
cprev = g.centroid(c,:);
gs{1} = g;
for i=1:sim_length
    filename = [fold,'\',strrep(name_for, 'sss', num2str(i))];
    load(filename,'g');
    for k=1:length(g.cells)-1
        if ~g.dead(k), g.centroid(k,:) = cellCOM(g,k); end
    end
    com = g.centroid(c,:);
    gs{i+1} = g;
    if norm(com-cprev) > jthresh && norm(com-cprev) < 6 , marked_frames = [marked_frames, i]; end
    cprev = com;
end

%% correct jumps
for i=1:length(marked_frames)
    mf = marked_frames(i);
    g1 = gs{mf}; % frame previous to the jump
    g2 = gs{mf+1}; % upload jumped frame
    % calculate the mean displacement of all the cells
    not_dead = (~g1.dead & ~g2.dead);
%     dx = mean(g2.centroid(not_dead,1) - g1.centroid(not_dead,1));
    dx = g2.centroid(c,1) - g1.centroid(c,1);
    % move every frame after the jump back with dx
    for j=mf:sim_length
        g = gs{j+1};
        g.verts(:,1) = g.verts(:,1) - dx;
        g = confineVerts(g);
        gs{j+1} = g;
    end
    disp(num2str(i)); % temp
end

%% make new movie
vid = VideoWriter('E:\Roie\vertex_code\no_jumps1.avi');
vid.Quality = 100; % 100 for best
vid.FrameRate = 12;
open(vid);
figure(10);
for t=0:sim_length
    g = gs{t+1};
    LatticePresentation(g,0, 10);
    lm = 3.7;
    xlim([-lm lm]);
    ylim([-lm lm]);
    frame = getframe;
    writeVideo(vid,frame);
end
close(vid);