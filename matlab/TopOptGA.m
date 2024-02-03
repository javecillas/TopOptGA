%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Controlling Out-of-Plane Buckling in Shear-Acting Structural Fuses
%%%%%% Through Topology Optimization
%%%%%% Javier A. Avecillas; Matthew R. Eatherton
%%%%%% Department of Civil and Environmental Engineering, Virginia Tech
%%%%%% Version 1.0 - Last update: 07/09/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Summary
% This routine performs the topology optimization of shear-acting
% structural fuses using genetic algorithms.
%   - The design domain has been discretized using 32 by 32 square elements
%   - Two extra rows of 1 by 32 elements are added to the top and bottom
%     boundary to facility the application of boundary conditions
%   - To guarantee the doubly-symmetry of the optimized topology, the
%     genetic algorithm only works with one quadrant of the actual topology
%   - The genetic operators are listed as follows: exponential ranking
%     selection, circular single-point crossover, constant volume mutation
%     and elitism condition.
%   - The units used in this routine are: inches, kips, and ksi

%% Clear workspace variables and close windows
clear all
close all
clc

%% General input
% Geometric dimensions of the domain (design domain + boundary elements)
width = 28.00;
height = 29.75;
thickness = 0.50;
% Number of finite elements used to discretize the domain
    % These numbers must be selected to guarantee that the domain is
    % discretized with SQUARE elements
    % Number of finite elements in the x direction
nx = 32;
    % Number of finite elements in the y direction - 2 extra rows for
    % boundary elements
ny = 34;
% Material properties - Steel A572 Gr 50
    % Young's modulus
E = 29000;
    % Yield strength
Fy = 50.0;
    % Poisson' ratio = 0.3 (For all calculations)
% Reference load for eigen-buckling analysis 
FEAP_ref_load = 1;
% Reference displacement for shear yield analysis 
FEAP_ref_disp = 0.595;

%% Genetic algorithm input
% Number of topologies in each generation
ini_pop = 100;
% Total number of generations
gnts = 400;
% Stop criteria
    % If 'stop_cri' >= 1 then run all 'gnts'
stop_cri = 1.1;
c_stop = 0;
% Volume fraction
    % The number of active elements in one quadrant is equal to 'n_ele_mat'
n_ele_mat = 102;
% Define upper-right quadrant information
    % The genetic algorithm only operates with one quadrant to enforce
    % doubly-symmetry of the optimized topology
nx_SD = 16;
ny_SD = 16;
% Connectivity parameter
    % Edges only - conn = 4
    % Edges and corners - conn = 8
conn = 4;
% Repair parameters
    % Invoked only if the largest component of the topology contains at
    % least 'min_ra' of the target volume fraction
min_ra = 0.90;
    % Maximum number of iterations
n_repair = 500;

%% Genetic algorithm operators
% Number of topologies considered as elite - must be an even number
n_eli = 4;
% Number of topologies considered for selection operation
n_sel = ini_pop-n_eli;
% Number of elements selected for crossover operation - 50% of total
xvr_ele = 128;
% Mutation operation
    % Probability of mutation
p_mut = 0.05;
    % Mutation area
m_area = 0.10;

%% Plot parameters
% Selected generations to be analized
gnts_vect = [1 2 5 10 20 50 100 200 300 400];
% Subplot configuration - #plots_x by #plots_y
plot_conf = [3,4];

%% Save generated information
% Saves all the generated topologies at each generation in a 3D array
shapes = zeros(nx_SD*ny_SD,ini_pop,gnts+1);
% Saves all the generated data at each generation in a 3D array
data = zeros(10,ini_pop,gnts+1);

%% Load initial population
% Important note: The txt file with the initial population will be deleted
% after the evaluation of the intial population
initial_popu = load('GA_IP.txt');

%% Evaluation of the intial population
for c_ip = 1:ini_pop
    % Run the repair algorithm
    [ R_Input, ~ ] = Repair_Input( initial_popu(:,c_ip), nx_SD, ny_SD, conn, min_ra, n_repair );
    % Construct the complete design domain of the topology
        % symm = 3 applies doubly-symmetry
    symm = 3;
    [ S_R_Input ] = Symmetry_Input( R_Input, nx_SD, ny_SD, symm );
    % Add top and bottom boundary elements
    Buck_B_S_R_Input = vertcat(ones(2*nx_SD,1),S_R_Input,ones(2*nx_SD,1));
    Yield_B_R_Input = vertcat(R_Input,ones(nx_SD,1));
    % Run the image-processing-based analysis to check load path
    [ data(4:9,c_ip,1) ] = Img_Processing( Buck_B_S_R_Input, width, height, nx, ny, nx_SD, ny_SD, conn );
    % Run FEAP FE analysis - if topology is admissible
    if data(5,c_ip,1) == 1 && data(6,c_ip,1) == 0
        % Shear yield FE analysis - Vy
        [ fetopo,fecoord ] = Generate_Mesh( width/2, height/2, nx/2, ny/2 );
        [ data(1,c_ip,1) ] = FEAP_Shear_Yielding( Yield_B_R_Input, width/2, height/2, thickness, nx/2, ny/2, E, Fy, FEAP_ref_disp/2, fetopo, fecoord );
        % Shear buckling FE analysis - Vb
        [ fetopo,fecoord ] = Generate_Mesh( width, height, nx, ny );
        [ data(2,c_ip,1) ] = FEAP_Shear_Buckling( Buck_B_S_R_Input, thickness, nx, ny, E, FEAP_ref_load, fetopo, fecoord );        
        % Compute radial distance sqrt(Vy^2+Vb^2)
        data(3,c_ip,1) = (data(1,c_ip,1)^2+data(2,c_ip,1)^2)^0.5;
        % Saves current topology
        shapes(:,c_ip,1) = R_Input(:,1);
    else
        % Saves current topology
        shapes(:,c_ip,1) = R_Input(:,1);        
    end
    % Delete all txt files
    delete *.txt
end
% Compute the objective function using Vy and Vb
[ data(10,:,1) ] = Objective_Function( data(:,:,1), 4*n_ele_mat, ini_pop );

%% Perform the genetic algorithm
for c_gnt = 1:400
    %% Print current generation
    c_gnt
    
    %% Run the selection operation
    Selection_OF = data(10,:,c_gnt);
    [ S_couples_idx ] = Selection_EXPRANK( Selection_OF, n_sel );
    
    %% Run the crossover operation
    child = zeros(nx_SD*ny_SD,n_sel);
    for c_xvr = 1:n_sel/2
        [ child(:,1+2*(c_xvr-1)), child(:,2+2*(c_xvr-1)) ] = Crossover_CSP( shapes(:,S_couples_idx(c_xvr,1),c_gnt), shapes(:,S_couples_idx(c_xvr,2),c_gnt), nx_SD, ny_SD, xvr_ele, 0, 0 );
    end
    
    %% Run the mutation operation
    for c_mut = 1:n_sel
        [ child(:,c_mut) ] = Mutation_EQVOL( child(:,c_mut), p_mut, m_area );
    end
    
    %% Run the elitism condition
    [ ~, idx_eli ] = sort(data(10,:,c_gnt));
    for c_eli = 1:n_eli
        child = horzcat(child,shapes(:,idx_eli(c_eli),c_gnt));
    end

    %% Evaluation of the current population
    for c_cp = 1:ini_pop
        % Run the repair algorithm
        [ R_Input, ~ ] = Repair_Input( child(:,c_cp), nx_SD, ny_SD, conn, min_ra, n_repair );
        % Construct the complete design domain of the topology
            % symm = 3 applies doubly-symmetry
        symm = 3;
        [ S_R_Input ] = Symmetry_Input( R_Input, nx_SD, ny_SD, symm );
        % Add top and bottom boundary elements
        Buck_B_S_R_Input = vertcat(ones(2*nx_SD,1),S_R_Input,ones(2*nx_SD,1));
        Yield_B_R_Input = vertcat(R_Input,ones(nx_SD,1));
        % Run the image-processing-based analysis to check load path
        [ data(4:9,c_cp,c_gnt+1) ] = Img_Processing( Buck_B_S_R_Input, width, height, nx, ny, nx_SD, ny_SD, conn );
        % Run FEAP analysis - if topology is admissible
        if data(5,c_cp,c_gnt+1) == 1 && data(6,c_cp,c_gnt+1) == 0
            % Shear yield FE analysis - Vy
            [ fetopo,fecoord ] = Generate_Mesh( width/2, height/2, nx/2, ny/2 );
            [ data(1,c_cp,c_gnt+1) ] = FEAP_Shear_Yielding( Yield_B_R_Input, width/2, height/2, thickness, nx/2, ny/2, E, Fy, FEAP_ref_disp/2, fetopo, fecoord );
            % Shear buckling FE analysis - Vb
            [ fetopo,fecoord ] = Generate_Mesh( width, height, nx, ny );
            [ data(2,c_cp,c_gnt+1) ] = FEAP_Shear_Buckling( Buck_B_S_R_Input, thickness, nx, ny, E, FEAP_ref_load, fetopo, fecoord );        
            % Compute radial distance sqrt(Vy^2+Vb^2)
            data(3,c_cp,c_gnt+1) = (data(1,c_cp,c_gnt+1)^2+data(2,c_cp,c_gnt+1)^2)^0.5;
            % Saves current topology
            shapes(:,c_cp,c_gnt+1) = R_Input(:,1);
        else
            % Saves current topology
            shapes(:,c_cp,c_gnt+1) = R_Input(:,1);
        end
    end
    % Compute the objective function using Vy and Vb
    [ data(10,:,c_gnt+1) ] = Objective_Function( data(:,:,c_gnt+1), 4*n_ele_mat, ini_pop );
    % Check stop criteria
    if min(data(10,:,c_gnt)) == min(data(10,:,c_gnt+1));
        c_stop = c_stop+1;
        if c_stop == round(stop_cri*gnts)
            break
        else
        end
    else
        c_stop = 0;
    end
    % Delete all txt files
    delete *.txt
end

%% Postprocessing
% Generate scatter plots of (Vb,Vb) at each 'gnts_vect'
[ f1, f2, f3, f4 ] = Plot_GA_Results( data, shapes, nx_SD, ny_SD, nx, ny, gnts_vect, plot_conf );
% Plot the fraction of admissible topologies at each 'gnts_vect'
[ f5 ] = Plot_Feasible( data );

%% Save the Shapes and Data as .mat files
save('GA_Data_400.mat','data')
save('GA_Shapes_400.mat','shapes')