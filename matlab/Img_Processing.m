%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Controlling Out-of-Plane Buckling in Shear-Acting Structural Fuses
%%%%%% Through Topology Optimization
%%%%%% Javier A. Avecillas; Matthew R. Eatherton
%%%%%% Department of Civil and Environmental Engineering, Virginia Tech
%%%%%% Version 1.0 - Last update: 07/09/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% IMAGE-PROCESSING-BASED ANALYSIS %%%%%%%%%%%%%%%%%%%%%
% 'Input':        2D binary array representing the repaired topology
%                 including the top and bottom boundary elements
% 'width':        Total width of the full domain
% 'height':       Total height of the full domain
% 'nx' and 'ny'   Number of elements in the x and y direction - Full domain
% 'nx_SD'         Number of elements in the x - Quadrant domain
% 'ny_SD'         Number of elements in the y - Quadrant domain
% 'conn'          Connectivity parameter
%                 Edges only - conn = 4
%                 Edges and corners - conn = 8

function [ img_prc ] = Img_Processing( Input, width, height, nx, ny, nx_SD, ny_SD, conn )
%% Reshape input vector into a 2D array
Input = reshape(Input,nx,ny)';

%% Pixel area
a_pix = (width/nx)*(height/ny);

%% Actual active elements
act_ae = (sum(Input(:))-4*nx_SD);

%% Check load path
% If l_p = 1 there is load path
% If l_p = 0 there is not load path
ele_lbl = bwlabel(Input,conn);
if ele_lbl(1,1) == ele_lbl(end,1)
    l_p = 1;
else
    l_p = 0;
end

%% Number of structural disconnected components
n_dcomp = max(ele_lbl(:))-l_p;

%% Area of structural disconnected components
if n_dcomp == 0;
    a_dcomp = [0 0];
else
    a_dcomp = zeros(n_dcomp,2);
    a_dcomp(1:end,1) = linspace(1+l_p,n_dcomp+l_p,n_dcomp);
    for cont_dcomp = 1:n_dcomp
        idx_dcomp = a_dcomp(cont_dcomp,1);
        a_dcomp(cont_dcomp,2) = sum(ele_lbl(:) == idx_dcomp)*a_pix;
    end
end
ta_dcomp = sum(a_dcomp(:,2));

%% Minimum area of disconnected components
mina_dcomp = min(a_dcomp(:,2));

%% Minimum area of internal holes
% Add external frame
holes = ones(size(Input,1)+2,size(Input,2)+2);
holes(2:end-1,2:end-1) = (1-Input);
% Label the holes present in the image
hole_lbl = bwlabel(holes,conn);
% Get the number of internal holes
n_holes  = max(hole_lbl(:))-1;
% Minimum area of internal holes
if n_holes == 0;
    a_holes = [0 0];
else
    a_holes = zeros(n_holes,2);
    a_holes(1:end,1) = linspace(2,n_holes+1,n_holes);
    for cont_h = 1:n_holes
        a_holes(cont_h,2) = sum(hole_lbl(:) == cont_h+1)*a_pix;
    end    
end
mina_ihol = min(a_holes(:,2));

%% Gather all the information
img_prc = [act_ae, l_p, n_dcomp, ta_dcomp, mina_dcomp, mina_ihol]';

end