% Pfad-Initialisierung für die Sammlung externer Matlab-Toolboxen

% Moritz Schappler, schappler@imes.uni-hannover.de, 2018-03
% (C) Institut für mechatronische Systeme, Universität Hannover

this_tb_path = fileparts( mfilename('fullpath') );
addpath(this_tb_path);

% Prüfe, ob die IMES-Robotik-Toolbox schon im Pfad ist.
% imes-robotics sollte erst nach matlab-ext initialisiert werden, damit
% Funktionen wie "eul2r" von dort genommen werden.
if ~isempty(which('robotics_toolbox_path_init.m'))
  warning('IMES-Robotik-Repo war schon im Pfad. Dadurch werden die Funktionen dort überlagert.');
end
addpath(fullfile(this_tb_path, 'miscellaneous'));
addpath(fullfile(this_tb_path, 'Advanced_Setpoints'));
addpath(fullfile(this_tb_path, 'export_fig'));
addpath(fullfile(this_tb_path, 'matGeom', 'geom2d'));
addpath(fullfile(this_tb_path, 'matGeom', 'geom3d'));
addpath(fullfile(this_tb_path, 'matGeom', 'meshes3d'));
addpath(fullfile(this_tb_path, 'matGeom', 'polygons2d'));
addpath(fullfile(this_tb_path, 'matGeom', 'utils'));
addpath(fullfile(this_tb_path, 'rvctools'));

run(fullfile(this_tb_path, 'rvctools', 'startup_rvc.m'));
% Ordner aus Pfad entfernen, die andere Skripte stören 
% (`codegen` überschattet die Matlab-interne Funktion im Pfad)
rmpath(fullfile(this_tb_path, 'rvctools', 'robot', 'demos'));
rmpath(fullfile(this_tb_path, 'rvctools', 'robot', 'examples'));
rmpath(fullfile(this_tb_path, 'rvctools', 'robot', 'mex'));
