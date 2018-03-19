% Pfad-Initialisierung für die Sammlung externer Matlab-Toolboxen

% Moritz Schappler, schappler@imes.uni-hannover.de, 2018-03
% (C) Institut für mechatronische Systeme, Universität Hannover

this_tb_path = fileparts( mfilename('fullpath') );
addpath(this_tb_path);

addpath(fullfile(this_tb_path, 'miscellaneous'));
addpath(fullfile(this_tb_path, 'Advanced_Setpoints'));
addpath(fullfile(this_tb_path, 'export_fig'));
addpath(fullfile(this_tb_path, 'geom2d', 'geom2d'));
addpath(fullfile(this_tb_path, 'geom2d', 'polygons2d'));
addpath(fullfile(this_tb_path, 'geom2d', 'polynomialCurves2d'));
addpath(fullfile(this_tb_path, 'geom2d', 'utils'));
addpath(fullfile(this_tb_path, 'geom3d'));
addpath(fullfile(this_tb_path, 'meshes3d'));
addpath(fullfile(this_tb_path, 'rvctools'));
run(fullfile(this_tb_path, 'rvctools', 'startup_rvc.m'));
