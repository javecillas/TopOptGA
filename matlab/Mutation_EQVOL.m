%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Controlling Out-of-Plane Buckling in Shear-Acting Structural Fuses
%%%%%% Through Topology Optimization
%%%%%% Javier A. Avecillas; Matthew R. Eatherton
%%%%%% Department of Civil and Environmental Engineering, Virginia Tech
%%%%%% Version 1.0 - Last update: 07/09/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% UNIFORM MUTATION OPERATION %%%%%%%%%%%%%%%%%%%%%%%
%%% This mutation operation preserves the volume of the topology
% 'Input':        Binary input vector
% 'p_mut':        Probability of mutation
% 'm_area':       Fraction of the active area to be mutated

function [ mut_Input ] = Mutation_EQVOL( Input, p_mut, m_area )
%% Step 1
% Get all the initial information
act_ele = sum(Input);
mut_ele = round(m_area*act_ele);
% Find indixes corresponding to active and void elements elements
[idx_act,~] = find(Input == 1);
[idx_void,~] = find(Input == 0);

%% Step 2
% Mutation procedure
r = rand;
if r <= p_mut && mut_ele > 0
    mut_Input = Input;
    % Change active elements into voids
    act_mut = datasample(idx_act,mut_ele,'Replace',false);
    mut_Input(act_mut,1) = 0;
    % Change voids into active elements
    void_mut = datasample(idx_void,mut_ele,'Replace',false);
    mut_Input(void_mut,1) = 1;    
else
    mut_Input = Input;
end

end