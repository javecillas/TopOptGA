%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Controlling Out-of-Plane Buckling in Shear-Acting Structural Fuses
%%%%%% Through Topology Optimization
%%%%%% Javier A. Avecillas; Matthew R. Eatherton
%%%%%% Department of Civil and Environmental Engineering, Virginia Tech
%%%%%% Version 1.0 - Last update: 07/09/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GENERATE MESH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'width'         Width of the domain
% 'height'        Height of the domain
% 'nx' and 'ny'   Number of elements in the x and y direction

function [ fetopo,fecoord ] = Generate_Mesh( width,height,nx,ny )
%% Set working directory and set path
hdir = pwd;
cd(hdir)
addpath(genpath(hdir));

%% Create new temp directory
tmpdir = pwd;
fid = fopen('tmpdir.log','a');
fprintf(fid,'Temporary directory: %s at %s\n\n',tmpdir,datestr(now));
fclose(fid);
cd(tmpdir);

%% Inputs for generating structured mesh
% Domain dimension
tDim  = 2;
% The x-coords of corners: [xmin xmax xmax xmin]
xcorn = [0.0 width width 0.0];
% The y-coords of corners: [xmin ymin ymax ymax] 
ycorn = [0.0 0.0 height height];
% Length in z-dir
meshpro(1,1) = 0.0;
% Number of elements in x-dir
meshpro(1,2) = nx;                 
% Number of elements in y-dir
meshpro(1,3) = ny;
% Number of elements in z-dir
meshpro(1,4) = 0;
if     ( tDim == 2 )
    [fetopo,fecoord] = MeshGenerate2D(xcorn,ycorn,meshpro);
elseif ( tDim == 3 )
    [fetopo,fecoord] = MeshGenerate3D(xcorn,ycorn,meshpro);
end

end

