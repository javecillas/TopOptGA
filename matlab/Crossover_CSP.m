%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Controlling Out-of-Plane Buckling in Shear-Acting Structural Fuses
%%%%%% Through Topology Optimization
%%%%%% Javier A. Avecillas; Matthew R. Eatherton
%%%%%% Department of Civil and Environmental Engineering, Virginia Tech
%%%%%% Version 1.0 - Last update: 07/09/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% CIRCULAR SINGLE POINT CROSSOVER OPERATION %%%%%%%%%%%%%%%%%
%%% This crossover operation preserves the volume of the topology
% 'parent_'       The parent topologies that will produce offspring
% 'nx' and 'ny'   Number of elements in the x and y direction
% 'alpha':        The number of [rows columns] to exchange in the crossover
%                 operator - set [0 0] for random choice
% 'p1' and 'p2':  Coordinates indicating the position where the cut
%                 begins - set [0 0] for random choice
% 'dir':          Directon of crossover =1 x-dir or =2 y-dir

function [ child_1, child_2 ] = Crossover_CSP( parent_1, parent_2, nx, ny, alpha, p_1, p_2 )
%% Direction of the crossover operation
dir = datasample([1 2],1);
if dir == 1;
elseif dir == 2;
    % Change direction of numbering for parent_1
    parent_1 = reshape(parent_1,nx,ny)';
    parent_1 = reshape(parent_1,nx*ny,1);
    % Change direction of numbering for parent_2
    parent_2 = reshape(parent_2,nx,ny)';
    parent_2 = reshape(parent_2,nx*ny,1);
end

%% Step 1
% Determine the number of elements
nele = length(parent_1);

%% Step 2
% Include the original position next to each allele
ori_p1 = horzcat(linspace(1,nele,nele)',parent_1);
ori_p2 = horzcat(linspace(1,nele,nele)',parent_2);

%% Step 3
% Determine initial 'p1' indicating the position where the cut begins in 'parent_1'
if isequal(p_1,0) == 1
    p_1 = randi([1 nele],1,1);
else   
end

%% Step 4
% Determine initial 'p2' indicating the position where the cut begins in 'parent_1'
if isequal(p_2,0) == 1
    p_2 = randi([1 nele],1,1);
else   
end

%% Step 5
% Shift parents' initial position to have p_1 and p_2 as first element
shi_p1 = circshift(ori_p1,-(p_1-1));
shi_p2 = circshift(ori_p2,-(p_2-1));

%% Step 6
% Check whether or not n1s_1 = n1s_2 and generate offspring
for cont_p1 = 0:nele-1
    cur_p1 = circshift(shi_p1,-cont_p1);
    cur_c1 = cur_p1;
    for cont_p2 = 0:nele-1
        cur_p2 = circshift(shi_p2,-cont_p2);
        cur_c2 = cur_p2;
        % Check for number of 1s
        n1s_1 = sum(cur_p1(1:alpha,2));
        n1s_2 = sum(cur_p2(1:alpha,2));
        eql = isequal(n1s_1,n1s_2);
        if eql == 1
            cur_c1(1:alpha,2) = cur_p2(1:alpha,2);
            cur_c2(1:alpha,2) = cur_p1(1:alpha,2);
            % Shift back vectors
            cur_c1 = circshift(cur_c1,cur_c1(1,1)-1);
            child_1 = cur_c1(:,2);
            cur_c2 = circshift(cur_c2,cur_c2(1,1)-1);
            child_2 = cur_c2(:,2);
            break
        elseif eql == 0 && cont_p1 == nele-1 && cont_p2 == nele-1
            cur_c1 = circshift(cur_c1,cur_c1(1,1)-1);
            child_1 = cur_c1(:,2);
            cur_c2 = circshift(cur_c2,cur_c2(1,1)-1);
            child_2 = cur_c2(:,2);
        end
    end
    if eql == 1
        break
    end
end

%% Post processing
% Change the order of the Input if direction Y was chosen
if dir == 1;
elseif dir == 2;
    % Change direction of numbering for parent_1
    child_1 = reshape(child_1,ny,nx);
    child_1 = reshape(child_1',nx*ny,1);
    % Change direction of numbering for parent_2
    child_2 = reshape(child_2,ny,nx);
    child_2 = reshape(child_2',nx*ny,1);
end

end