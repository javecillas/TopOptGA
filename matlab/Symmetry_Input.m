%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Controlling Out-of-Plane Buckling in Shear-Acting Structural Fuses
%%%%%% Through Topology Optimization
%%%%%% Javier A. Avecillas; Matthew R. Eatherton
%%%%%% Department of Civil and Environmental Engineering, Virginia Tech
%%%%%% Version 1.0 - Last update: 07/09/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% SYMMETRY CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This routine resembles the actual topology based on an binary input
%%% vector that represents only one quadrant of the topology
% 'input'         Binary input vector that represents one quadrant of the
%                 topology
% 'nx' and 'ny'   Number of elements in the x and y direction - Quadrant
% 'symm'          Symmetry condition
%                 symm = 1 -> Single symmetric about XX
%                 symm = 2 -> Single symmetric about YY
%                 symm = 3 -> Double symmetric                

function [ S_Input ] = Symmetry_Input( Input, nx, ny, symm )
%% Re-arrange binary input vector into a binary matrix
Input = flipud(reshape(Input,[nx,ny])');

%% Applyng symmetry conditions
% Doubly-symmetry always enforced with symm = 3
if symm == 1;
    S_Input = vertcat(Input,flipud(Input));
elseif symm == 2;
    S_Input = horzcat(fliplr(Input),Input);
    elseif symm == 3;
        S_Input = vertcat(Input,flipud(Input));
        S_Input = horzcat(fliplr(S_Input),S_Input);
end
%% Reshape binary array into a binary vector
S_Input = reshape(S_Input',size(S_Input,1)*size(S_Input,2),1);

end

