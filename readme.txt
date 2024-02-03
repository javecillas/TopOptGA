%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Controlling Out-of-Plane Buckling in Shear-Acting Structural Fuses
%%%%%% Through Topology Optimization
%%%%%% Javier A. Avecillas; Matthew R. Eatherton
%%%%%% Department of Civil and Environmental Engineering, Virginia Tech
%%%%%% Version 1.0 - Last update: 07/09/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Summary
This routine performs the topology optimization of shear-acting
structural fuses using genetic algorithms
	- The units used in this routine are: inches, kips, and ksi
	- The design domain has been discretized using 32 by 32 square elements
	- Two extra rows of 1 by 32 elements are added to the top and bottom
	boundary to facility the application of boundary conditions, but they
	are not part of the optimization routine
	- To guarantee the doubly-symmetry of the optimized topology, the
    genetic algorithm only operates with one quadrant of the actual topology
	- The genetic operators are listed as follows: exponential ranking
    selection, circular single-point crossover, constant volume mutation
    and elitism condition.

%% Example initial population
	- The initial population consists of binary input vectors representing
	one quadrant of the actual topology, excluding the top and bottom
	boundary elements
	- An example of an initial population with a volume fraction of 40%
	is provided in '/matlab/GA_IP.txt'
	- The 'GA_IP.txt' file must be in the same directory as the Matlab
	routines

%% Finite element software
	- This topology optimization routine uses FEAP as the finite element
	software, http://projects.ce.berkeley.edu/feap/
	- The FEAP executable must be in the same directory as the Matlab
	routines

%% Steps
The steps for using this routine are:
	1. Open the File 'TopOptGA.m'
	2. Input the required information including:
		2.1 General input, line 27 to 47
			- Geometric dimensions of the domain, including boundary elements
			- Number of finite elements used to discretize the domain
			- Material properties
			- Reference load and maximum displacement
		2.2 Genetic algorithm input, line 51 to 77
			- Number of topologies
			- Number of generations
			- Volume fraction
			- Connectivity condition
			- Repair procedure
		2.3 Genetic algorithm operators input, line 79 to 90
			- Number of elite topologies, must be an even number
			- Number of elements selected for crossover, 50% of total
			- Probability of mutation and mutation area
		2.4 Plot parameters input, line 92 to 96
			- Selected generations for post-processing
			- Plot configuration
	3. Load the initial population, 'GA_IP.txt'
	4. Evaluate the initial population
		4.1 Run repair algorithm
		4.2 Construct the complete topology
		4.3 Add top and bottom boundary elements
		4.4 Run the image-processing-based analysis to check load path
		4.5 Run shear yield FE analysis, get Vy
		4.6 Run shear buckling FE analysis, get Vb
		4.7 Compute and save the objective function value
		4.8 Save the current topology
	5. Perform the genetic algorithm
		5.1 Run the selection operation
		5.2 Run the crossover operation
		5.3 Run the mutation operation
		5.4 Run the elitism condition
		5.5 Evaluate the current population
			5.5.1 Run repair algorithm
			5.5.2 Construct the complete topology
			5.5.3 Add top and bottom boundary elements
			5.5.4 Run the image-processing-based analysis to check loadpath
			5.5.5 Run shear yield FE analysis, get Vy
			5.5.6 Run shear buckling FE analysis, get Vb
			5.5.7 Compute and save the objective function value
			5.5.8 Save the current topology
		5.6 Check stop condition
	6. Post-processing
		6.1 Plot Vb vs Vy and other plots using the generated information
		6.2 Plot the fraction of admissible topologies
	7. Save the files as '.mat' files
		7.1 Save the date as 'GA_Data_400.mat'
			- This file contains all the numerical information for each 
			topology at each generation
		7.2 Save the generated topologies as 'GA_Shapes_400.mat'
			- This file contains all the topologies at each generation
			- The topologies are binary vectors representing one quadrant
			of the actual topology
			