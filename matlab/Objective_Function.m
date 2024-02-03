%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Controlling Out-of-Plane Buckling in Shear-Acting Structural Fuses
%%%%%% Through Topology Optimization
%%%%%% Javier A. Avecillas; Matthew R. Eatherton
%%%%%% Department of Civil and Environmental Engineering, Virginia Tech
%%%%%% Version 1.0 - Last update: 07/09/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% OBJECTIVE FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'data'            Information about each topology at each generation
%                   data(1,i) -> Vy
%                   data(2,i) -> Vb
%                   data(3,i) -> sqrt(Vy^2+Vb^2)
%                   data(4,i) -> Number of active components
%                   data(5,i) -> Load path condition 1-yes 0-no
%                   data(6,i) -> Number of disconnected components
%                   data(7,i) -> Area of disconnected components
%                   data(8,i) -> Min. area of disconnected components
%                   data(9,i) -> Min. area of internal holes
%                   data(10,i) -> Objective function value
% 'n_ele_mat'       The number of active elements in one quadrant
%                   (without boundary elements)
% 'n_ele_mat'       The number of topologies in each generation

function [ data_EOF ] = Objective_Function( data, n_ele_mat, n_ind )
%% Assign weights for the penalty term - Normal distribution shape
% Maximum value of the penalty term when Vy/Vb = k
a_g = 1.00;
% Target value for Vy/Vb
b_g = 0.40;
% Width - We recomend values ~ 0.1 
% If c_g < 0.1 the penalty will be more severe and may cuase troubles
c_g = 0.10;

%% Assign values for equal area penalty
% If we preserve the volume fraction of all the topologies during the
% genetic algorithm, this penalty term is useless
AP1 = 1.0;

%% Assign values for connectivity penalty (not used in this study)
% Penalty for the total number of disconnected components
CP1 = 0.0;
% Penalty for the total area of disconnected components
CP2 = 0.0;
% Penalty for the minimum area of disconnected components
CP3 = 0.0;

%% Assign values for internal hole's penalty  (not used in this study)
% Set the minimum admisible area for internal holes
AdmA_Int_Hole = 1.0;
% Penalty for the minium area of internal holes
HP1 = 0.0;

%% Find indixes corresponding to OK topologies w/o disconnected components
[~,idx_OKstr] = find(data(2,:) ~= 0 & data(5,:) == 1 & data(6,:) == 0);

%% Find remaining indixes
idx_NOstr = linspace(1,n_ind,n_ind);
idx_NOstr(idx_OKstr) = [];

%% Objective function for OK topologies w/o disconnected components
data_OF = zeros(1,n_ind);
data_EOF = zeros(1,n_ind);
if isempty(idx_OKstr) == 1;
else
    for c_ind = 1:length(idx_OKstr)
        % Shear yield load, Vb
        Vy = data(1,idx_OKstr(c_ind));
        % Shear Buckling load, Vb
        Vb = data(2,idx_OKstr(c_ind));
        % Position vector, r
        r = data(3,idx_OKstr(c_ind));
        % Total active elements
        TA_Act_Ele = data(4,idx_OKstr(c_ind));
        % Load path | 1 - Yes 0 - No
        LP = data(5,idx_OKstr(c_ind));        
        % Number of structural disconnected components
        N_Str_DComp = data(6,idx_OKstr(c_ind));
        % Total area of structural disconnected components
        TA_Str_DComp = data(7,idx_OKstr(c_ind));
        % Minimum area of structural disconnected components
        MinA_Str_DComp = data(8,idx_OKstr(c_ind));
        % Minimum area of internal holes
        MinA_Int_Hole = data(9,idx_OKstr(c_ind));
        
        % Objective function
        OF = 1/(a_g*exp(-(Vy/Vb-b_g)^2/(2*c_g^2))*r);
        % Storage the value of the objective function
        data_OF(1,idx_OKstr(c_ind)) = OF;
        
        % Enhanced objective function
        % Penalties
        AP = AP1*(exp(abs((TA_Act_Ele-n_ele_mat)/n_ele_mat)));
        CP = CP1*N_Str_DComp+CP2*TA_Str_DComp+CP3*MinA_Str_DComp;
        HP = HP1*(AdmA_Int_Hole-MinA_Int_Hole);
        % Calculation of EOF
        EOF = AP*OF+CP+HP;
        % Storage the value of the Enhanced Objective Function
        data_EOF(1,idx_OKstr(c_ind)) = EOF;
    end
end

%% Compute the objective function for the remaining of the population
% Get the worst objective function from feasible individuals
if isempty(idx_OKstr) == 1;
    worst_OF = 1;
else
    worst_OF = 1.01*max(data_OF(:));
end
% Compute the objective function for non feasible individuals
if isempty(idx_NOstr) == 1;
else
    for c_ind = 1:length(idx_NOstr)
        % Total active elements
        TA_Act_Ele = data(4,idx_NOstr(c_ind));
        % Number of structural disconnected domponents
        N_Str_DComp = data(6,idx_NOstr(c_ind));
        % Total area of structural disconnected domponents
        TA_Str_DComp = data(7,idx_NOstr(c_ind));
        % Minimum area of structural disconnected domponents
        MinA_Str_DComp = data(8,idx_NOstr(c_ind));
        % Minimum area of internal holes
        MinA_Int_Hole = data(9,idx_NOstr(c_ind));
        
        % Enhanced objective function
        % Penalties
        AP = AP1*(exp(abs((TA_Act_Ele-n_ele_mat)/n_ele_mat)));
        CP = CP1*N_Str_DComp+CP2*TA_Str_DComp+CP3*MinA_Str_DComp;
        HP = HP1*(AdmA_Int_Hole-MinA_Int_Hole);
        % Calculation of EOF
        EOF = AP*worst_OF+CP+HP;
        % Storage the value of the Enhanced Objective Function
        data_EOF(1,idx_NOstr(c_ind)) = EOF;        
    end
end

end