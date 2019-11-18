% run over all the lattices in a specific folder

fold_name = 'simulations-191118';
lat_files = dir([fold_name, '/lat_step(20)*.mat']);
nl = length(lat_files);
gs = cell(nl);
% load all lattices
for i=1:nl
    load([lat_files(i).folder, '/', lat_files(i).name], 'g');
    gs{i} = g;
end

proper = ones(nl,1); %logical array for proper lattices
% check for lattices with breakage in pillar row
% if we find a PC with only one PC neighbor, the row is broken.
for i=1:nl
    g = gs{i};
    PCs = find(g.populations == 4);
    for ci = 1:length(PCs)
        c = PCs(ci);
        if sum( g.populations(g.bonds(g.cells{c+1},4)) == 4 ) <= 1 
            proper(i) = 0;
            break;
        end
    end
end

%% mark best lattices
best_lat = zeros(nl,1);
for i=1:nl
    g = gs{i};
    if proper(i)
        LatticePresentation(g,0);
        [~, ~, bot] = ginput2(1);
        if (bot == 115) || (bot==1) % left MC or 's'
            best_lat(i) = 1;
        end
    else
        continue;
    end
    clf
end
    
%% mark best of the best lattices
bb_lat = zeros(nl,1);
for i=1:nl
    g = gs{i};
    if best_lat(i)
        LatticePresentation(g,0);
        [~, ~, bot] = ginput2(1);
        if (bot == 115) || (bot==1) % left MC or 's'
            bb_lat(i) = 1;
        end
    else
        continue;
    end
    clf
end

%%
 i=0;
 bot = 0;
 while bot~=27
     load(['simulations-191118/lat_step(', num2str(i), ')_comp(6)_tens(0.3)_cont(0.15)_pp(0.5)_sh(0.4).mat'],'g');
     LatticePresentation(g,0)
     set(gca,'Visible','off');
     [~, ~, bot] = ginput2(1);
     switch bot
         case 28
             i = mod(i-1,21);
         case 29
             i = mod(i+1, 21);
     end
     clf;
 end
