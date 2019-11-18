
% runStage1;
% waitbar(0.4,'Stage 2');
% runStage2;
% waitbar(0.8,'Statistics');
% produce_stats('stage1_2019-05-22',50,400,1);
% produce_stats('stage2_2019-05-22',50,100,2);
% plot_statistics;


% 'gl-alpha' - the global weight of the area term
% 'gl-gamma' - the global weight of the tension term
% 'gl-Gamma' - the global weight of the contractility term
% 'gl-eta' - the global weight of the shear term
% 'gl-zeta' - the global weight of the compression term
% 'gl-sigma' - the global weight of the repulsion term
% 'kappa' - the exponent of the repulsion term
% 'HC-Gamma' - the contractility of HCs
% 'HC-alpha' - the compressibility of HCs
% 'HCSCt' - the tension in HC:SC bonds

g.paras = [12 ; 0.3; 0.15; 0; 1; 1; 0.5];
sim_num = 20;
wb = waitbar(0,'Parameters run');
runStage1('2019-05-22', sim_num, 400, 'gl-alpha', (0.5:0.25:1.5)*12);
waitbar(0.12);
runStage1('2019-05-22', sim_num, 400, 'gl-gamma', (0.5:0.25:1.5)*0.3);
waitbar(0.24);
runStage1('2019-05-22', sim_num, 400, 'gl-Gama', (0.5:0.25:1.5)*0.15);
waitbar(0.36);
runStage1('2019-05-22', sim_num, 400, 'gl-eta', (0.5:0.25:1.5)*0.5);
waitbar(0.5);
runStage1('2019-05-22', sim_num, 400, 'gl-zeta', (0.5:0.25:1.5)*1);
waitbar(0.62);
runStage1('2019-05-22', sim_num, 400, 'gl-sigma', (0.5:0.25:1.5)*1);
waitbar(0.75);
runStage1('2019-05-22', sim_num, 400, 'HC-Gamma', (0.5:0.25:1.5)*2);
waitbar(0.87);
runStage1('2019-05-22', sim_num, 400, 'HC-alpha', (0.5:0.25:1.5)*2);
waitbar(1);


produce_stats_params('stage1_2019-05-22_gl-gamma', 20, 400, 'gl-gamma', (0.5:0.25:1.5)*0.3, 'stage1_2019-05-22_params-stats');
produce_stats_params('stage1_2019-05-22_gl-Gama', 20, 400, 'gl-Gama', (0.5:0.25:1.5)*0.15, 'stage1_2019-05-22_params-stats');
produce_stats_params('stage1_2019-05-22_gl-eta', 20, 400, 'gl-eta', (0.5:0.25:1.5)*0.5, 'stage1_2019-05-22_params-stats');
produce_stats_params('stage1_2019-05-22_gl-zeta', 20, 400, 'gl-zeta', (0.5:0.25:1.5)*1, 'stage1_2019-05-22_params-stats');
produce_stats_params('stage1_2019-05-22_gl-sigma', 20, 400, 'gl-sigma', (0.5:0.25:1.5)*1, 'stage1_2019-05-22_params-stats');
produce_stats_params('stage1_2019-05-22_HC-Gamma', 20, 400, 'HC-Gamma', (0.5:0.25:1.5)*2, 'stage1_2019-05-22_params-stats');
produce_stats_params('stage1_2019-05-22_HC-alpha', 20, 400, 'HC-alpha', (0.5:0.25:1.5)*2, 'stage1_2019-05-22_params-stats');
