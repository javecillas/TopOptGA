%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Controlling Out-of-Plane Buckling in Shear-Acting Structural Fuses
%%%%%% Through Topology Optimization
%%%%%% Javier A. Avecillas; Matthew R. Eatherton
%%%%%% Department of Civil and Environmental Engineering, Virginia Tech
%%%%%% Version 1.0 - Last update: 07/09/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% MESH GENERATION 2D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'xcorn'         The x corner points of the domain
% 'ycorn'         The y corner points of the domain
% 'meshpro'       (1,1) Length in z-dir
%                 (1,2) Number of elements in x-dir
%                 (1,3) Number of elements in y-dir
%                 (1,4) Number of elements in z-dir

function [fetopo,fecoord] = MeshGenerate2D(xcorn,ycorn,meshpro)
%% Get generation parameters
mshx = meshpro(1,2);
mshy = meshpro(1,3);

%% Mesh error control
if mshx < 1 || mshy < 1;
    display('mesh number should be > 1');
    error ('invalid mesh number or mesh number < 1');
end

%% Get mesh division
mshpsi = linspace(-1,1,mshx+1);
msheta = linspace(-1,1,mshy+1);

%% Get the coordinates of corner points
crpnt = [ xcorn(1,1)   ycorn(1,1)  ;...
          xcorn(1,2)   ycorn(1,2)  ;...
          xcorn(1,3)   ycorn(1,3)  ;...
          xcorn(1,4)   ycorn(1,4) ];

%% Get grid
[XI,ETA] = meshgrid(mshpsi,msheta);

%% Get the expanded shape function values
N1 = 1/4*(1-XI).*(1-ETA);
N2 = 1/4*(1+XI).*(1-ETA);
N3 = 1/4*(1+XI).*(1+ETA);
N4 = 1/4*(1-XI).*(1+ETA);

%% Get the coordinate of meshgrid
X = N1.*crpnt(1,1) + N2.*crpnt(2,1) + N3.*crpnt(3,1) + N4.*crpnt(4,1);
Y = N1.*crpnt(1,2) + N2.*crpnt(2,2) + N3.*crpnt(3,2) + N4.*crpnt(4,2);

%% Rearrange coordinate
X = (reshape(X',[],1));
Y = (reshape(Y',[],1));

%% Get the x,y,z coordinates in local and global
fecoord = [roundn(X,-14) roundn(Y,-14)];

%% Create connectivity matrix
no = 1;
for i = 1:mshy
    for j = 1:mshx
        if i < mshy + 1 && j < mshx + 1;  
            % Underside nodes
            fetopo(no,1) = (mshx+1)*(i-1)+j  ;
            fetopo(no,2) = (mshx+1)*(i-1)+j+1;
            fetopo(no,3) = (mshx+1)*(i)+j+1  ;
            fetopo(no,4) = (mshx+1)*(i)+j    ;
            
        end
        no = no + 1;
    end
end

end