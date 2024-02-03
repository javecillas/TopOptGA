%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Controlling Out-of-Plane Buckling in Shear-Acting Structural Fuses
%%%%%% Through Topology Optimization
%%%%%% Javier A. Avecillas; Matthew R. Eatherton
%%%%%% Department of Civil and Environmental Engineering, Virginia Tech
%%%%%% Version 1.0 - Last update: 07/09/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% EXPONENTIAL RANKING SELECTION %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% FOR A MINIMIZATION PROBLEM %%%%%%%%%%%%%%%%%%%%%%%%%
% 'input_OF'      Is a row vector containing the OF values of the
%                 topologies
% 'nsel'          Number of individuals to be selected for mating
% 'n_eta/nsel'    Probability of the worst candidate to be selected  
% 'couples_idx'   Indexes of randomly generated couples
% Exponent "c" (0 1) // Line 34 // Important

function [ couples_idx ] = Selection_EXPRANK( input_OF, nsel )
%% Step 1
% Identify the populaiton size
[~,N] = size(input_OF);

%% Step 2
% Reordering the Population OF in a descending order
input_OF = vertcat(linspace(1,N,N),input_OF);
[~,idx] = sort(input_OF(2,:),'descend');
sort_input_OF = input_OF(:,idx);

%% Step 3
% Create ranking positions
rank = linspace(1,N,N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Base of the exponent - "c" coefficient in the range (0 1)
% The closer to 1, the lower the exponentiality of the procedure
c = 0.80;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign a exponential distribution of probabilities to each position
prob = ((c-1)/(c^N-1))*c.^(N-rank);
% Normalize probabilities
norm_prob = prob/sum(prob);
% Cumulative probabilities
cum_prob = cumsum(norm_prob);
% Concatenate original positions and comulative probabilities
pot_parent = [sort_input_OF(1,:); cum_prob];

%% Step 4
% Use stochastic universal sampling to pick the corresponding parents
% Distance between pointers
p_dist = pot_parent(2,N)/nsel;
% Random start point
a = 0;
b = p_dist;
start = (b-a).*rand(1,1)+a;
% Calculate trials positions
trials = linspace(0,nsel-1,nsel);
trials = start+p_dist*trials;
% Find positions in the original numbering corresponding to the trials vals
Mf = repmat(pot_parent(2,:)',1,nsel);
Mt = repmat(trials,N,1);
[parent_idx, ~] = find(Mt<Mf & [zeros(1, nsel); Mf(1:N-1, :)]<=Mt);

%% Step 4
% Extract parents positions based in original numbering
parent_ori_idx  = pot_parent(1,parent_idx);
% Shuffle parents indexes
parent_ori_idx = parent_ori_idx(randperm(length(parent_ori_idx)));
% Create vector of couples
couples_idx = reshape(parent_ori_idx,[2,nsel/2])';

end