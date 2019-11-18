global timer

ext = '';
sf = '25'; %last file timestamp
datestr = '231218';
wb = waitbar(0,'winter is coming');
for i=1:50
    try
        load(['stage1_complete-',datestr,'\lat',ext,'(',num2str(i) ,')_step(',sf,')_diff.mat'],'g');
%         g = stage2differentiate(g);
        
        g.area_feedback = 0;
        g.fa = [1 1 1 1.2 1];
        g = redistributeAreas(g);
        
        timer = 0;
        % parameters: area, tension, perimeter, mechanical feedback, HCs vdv, PCs
        % pull, axial flow of SCs
        g.paras = [4 ; 0.4; 0.1; 0; 1; 0; 0];
        T1prob = 1;
        T1eps = 0.1;
        for t =1:31
            g = findTransitions(g,T1eps,T1prob,0.02); %0.05, 0.02, 0.02
            g = updateParameters(g);
            g = relaxLattice(g,5); %50
            
            g = killSmallCells(g,0.30);
            g = updateParameters(g);
            g = relaxLattice(g,5); %200
            
            disp(num2str(t))
            
            if mod(t-1,3) == 0
                save(['stage2_complete/lat',ext,'(', num2str(i),')_step(', num2str((t-1)/3) ,').mat'],'g');
            end
        end
    catch
        disp(['error in: ', num2str(i)]);
    end
    waitbar(i/50);
end
close(wb)

