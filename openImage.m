% selecting and loading the cia (cells image analysis) 
db_path = 'Q:\Users\Roie\database\';
filename = 'base 2_e17.5_28.1.mat';
% filename = 'base 3_e17.5_28.1.mat';
% filename = uigetfile(db_path);
load([db_path filename]);

imshow(a.Im);