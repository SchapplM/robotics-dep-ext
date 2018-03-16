% Gesamttest für die Matlab-Geometrie-Toolbox
% 
% Führt alle verfügbaren Modultests aus um die Funktionalität
% sicherzustellen

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-03
% (C) Institut für mechatronische Systeme, Universität Hannover

%% geom2d
clc;clear;close all;
this_repo_path = fullfile(fileparts(which('matlab_ext_test_repo.m')));
addpath(fullfile(this_repo_path,'examples_tests', 'geom2d','triangle'));
compTriangle
demoNapoleon
triangleDemo

%% geom3d
clc;clear;close all;
this_repo_path = fullfile(fileparts(which('matlab_ext_test_repo.m')));
addpath(fullfile(this_repo_path, 'examples_tests', 'geom3d-demos'));
demoDrawLine3d
demoDrawPlane3d
demoDrawTubularMesh
demoGeom3d
demoInertiaEllipsoid
demoRevolutionSurface
drawSoccerBall

%% meshes3d
clc;clear;close all;
this_repo_path = fullfile(fileparts(which('matlab_ext_test_repo.m')));
addpath(fullfile(this_repo_path, 'examples_tests', 'meshes3d-demos'));

% createTrefoilKnot
demoClipMeshVertices
demoConcatenateMeshes
demoCutMeshByPlane
demoPolyhedra
demoRemoveMeshFaces
demoRemoveMeshVertices
demoSplitMesh
demoTriangulateFaces
demoVoronoiCell

%% Ende
clc;clear;close all;
fprintf('Alle Testfunktionen dieses Repos ausgeführt\n');