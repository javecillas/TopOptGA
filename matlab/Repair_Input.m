%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Controlling Out-of-Plane Buckling in Shear-Acting Structural Fuses
%%%%%% Through Topology Optimization
%%%%%% Javier A. Avecillas; Matthew R. Eatherton
%%%%%% Department of Civil and Environmental Engineering, Virginia Tech
%%%%%% Version 1.0 - Last update: 07/09/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% TOPOLOGY REPAIRING %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'Input':        2D binary array representing the topology to be repaired
%                 withouth the top and bottom boundary elements
% 'nx' and 'ny'   Number of elements in the x and y direction
% 'conn'          Connectivity parameter
%                 Edges only - conn = 4
%                 Edges and corners - conn = 8
% 'min_a'         Fraction of the number of the total active elements that
%                 triggers the repair algorithm
% 'n_rep'         Maximum number of iterations

function [ R_Input, check_n_comp ] = Repair_Input( Input, nx, ny, conn, min_a, nrep )
for cont_rep = 1:nrep
    %% Re-arrange binary input vector into a binary matrix
    Input = (reshape(Input,[nx,ny]))';
    ele_num = linspace(1,nx*ny,nx*ny);
    ele_num  = reshape(ele_num,[nx,ny])';

    %% Identify the connected components and the dominant componet
    CC = bwconncomp(Input',conn);
    L = labelmatrix(CC)';
    n_comp = max(max(L));
    if n_comp == 1
        R_Input = Input;
        % Reshape binary array into a binary vector
        R_Input = reshape(R_Input',nx*ny,1);
        % Recheck number of subdomains
        Input = (reshape(R_Input,[nx,ny]))';
        check_CC = bwconncomp(Input',conn);
        check_L = labelmatrix(check_CC)';
        check_n_comp = max(check_L(:));
        % Break the iterative repairment
        return
    else
        C  = cell(n_comp,1);
        for i = 1:n_comp
            C(i,1) = CC.PixelIdxList(1,i);
        end
        num_ele = cellfun(@numel,CC.PixelIdxList);
        [~,idx_max] = max(num_ele);
        C_domin = CC.PixelIdxList{idx_max};
        % Check whether the dominant SUB domain meets the min area requirment
        if length(C_domin) < min_a*sum(num_ele)
            R_Input = Input;
            % Reshape binary array into a binary vector
            R_Input = reshape(R_Input',nx*ny,1);
            % Recheck number of subdomains
            Input = (reshape(R_Input,[nx,ny]))';
            check_CC = bwconncomp(Input',conn);
            check_L = labelmatrix(check_CC)';
            check_n_comp = max(check_L(:));
            % Break the iterative repairment
            return
        else
        end
        
        C_nondomin = CC.PixelIdxList;
        C_nondomin(idx_max) = [];

        %% Minimum distance between dominant component and remaining components
        D = 10^5*ones(n_comp-1,1);
        for i = 1:n_comp-1
            for j = 1:size(C_nondomin{i});
                for k = 1:size(C_domin,1);
                    Actual_C_nondomin = C_nondomin{i};
                    [x1,y1] = find(ele_num == Actual_C_nondomin(j,1));
                    [x2,y2] = find(ele_num == C_domin(k,1));
                    % Distance to the CG of each element
                    Dist = sqrt((x2-x1)^2+(y2-y1)^2);
                    if Dist < D(i,1)
                        D(i,1) = Dist;
                    else
                    end
                end
            end
        end

        %% Number of neighboring elements with actual material
        mask_NV = [];
        if conn == 4;
            mask_NV = [0 1 0; 1 0 1; 0 1 0];
        end
        if conn == 8;
            mask_NV = [1 1 1; 1 0 1; 1 1 1];
        end
        NV_Array = conv2(Input,mask_NV,'same');
        for i = 1:nx*ny
            if Input(i) == 1;
                NV_Array(i) = 0;
            else
            end
        end

        %% The sum of the weight of the neighbors
        NC_Array = zeros(ny,nx);
        for i = 1:nx*ny
            if Input(i) == 1;
            else
                sz = size(Input);
                % Row, Col index of the element
                [r,c] = ind2sub(sz,i);
                % Calculate of neighbors
                mask_NC = [];
                if conn == 4;
                    mask_NC(1:4,1:2) = [r+[0;-1;1;0] c+[-1;0;0;+1]];
                end
                if conn ==8;
                    mask_NC(1:8,1:2) = [r+[-1;0;1;-1;1;-1;0;1] c+[-1;-1;-1;0;0;1;1;1]];
                end
                % Select only those elements in the range
                mask_NC = mask_NC(all(mask_NC,2) & mask_NC(:,1)<=sz(1) & mask_NC(:,2)<=sz(2),:);
                % Convert to subindex notation to linear position
                idx = (mask_NC(:,2)-1)*sz(1) + mask_NC(:,1);
                NC_SubD = zeros(size(idx,1),1);
                for j = 1:size(idx,1)
                    NC_SubD(j,1) = L(idx(j));
                end
                % Delete group 0
                NC_SubD(NC_SubD == 0) = [];
                % Delete repeated sub domains
                NC_SubD = unique(NC_SubD);
                if isempty(NC_SubD) == 1
                    NC_Array(i) = 0;
                else
                    for j = 1:size(NC_SubD)
                        NC_Array(i) = NC_Array(i)+num_ele(NC_SubD(j));
                    end
                end
            end
        end

        %% Summation of NV and NC
        V_Array = NV_Array+NC_Array;

        %% Displace the most distant elements to the elements with greater V
        V_max = max(V_Array(:));
        [r_ava,c_ava] = find(V_Array == V_max);
        idx_ava = horzcat(r_ava,c_ava);
        n_ava_ele = size(idx_ava,1);

        farthest_SD = max(D(:));
        [r_farSD,~] = find(D == farthest_SD);
        mov_ele_list = [];
        for i = 1:size(r_farSD,1)
            mov_ele_list = vertcat(mov_ele_list,C_nondomin{r_farSD(i)});
        end
        n_mov_ele = size(mov_ele_list,1);
        idx_mov = [];
        for i = 1:n_mov_ele
            [r_mov,c_mov] = find(ele_num == mov_ele_list(i));
            idx_mov = vertcat(idx_mov,[r_mov,c_mov]);
        end
        if n_ava_ele >= n_mov_ele
            Rnd_r = randperm(size(idx_ava,1),n_mov_ele);
            idx_ava = idx_ava(Rnd_r,:);
        else
            Rnd_r = randperm(size(idx_mov,1),n_ava_ele);
            idx_mov = idx_mov(Rnd_r,:);
        end
        % Changing void and full elements
        R_Input = Input;
        for i = 1:min(n_mov_ele,n_ava_ele)
            % Changing 1s to 0s
            R_Input(idx_mov(i,1),idx_mov(i,2)) = 0;
            % Changing 0s to 1s
            R_Input(idx_ava(i,1),idx_ava(i,2)) = 1;
        end
        %% End of iterative process
        if cont_rep == nrep;
            % Reshape binary array into a binary vector for end of iteration
            R_Input = reshape(R_Input',nx*ny,1);
            % Recheck number of subdomains
            Input = (reshape(R_Input,[nx,ny]))';
            check_CC = bwconncomp(Input',conn);
            check_L = labelmatrix(check_CC)';
            check_n_comp = max(check_L(:));
        else
            % Reshape binary array into a binary vector for next iteration
            Input = reshape(R_Input',nx*ny,1);
        end        
    end
end

end