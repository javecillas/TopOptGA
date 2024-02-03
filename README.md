# **TopOptGA**

Avecillas, Javier A., and Matthew R. Eatherton. "Controlling out-of-plane
buckling in shear-acting structural fuses through topology optimization."
Journal of Structural Engineering 146.7 (2020): 04020132.

## Summary
This routine performs the topology optimization of shear-acting
structural fuses using genetic algorithms

- The units used in this routine are: inches, kips, and ksi.
- The design domain has been discretized using 32 by 32 square elements.
- Two extra rows of 1 by 32 elements are added to the top and bottom
boundary to facility the application of boundary conditions, but they
are not part of the optimization routine.
- To guarantee the doubly-symmetry of the optimized topology, the
genetic algorithm only operates with one quadrant of the actual topology.
- The genetic operators are listed as follows: exponential ranking
selection, circular single-point crossover, constant volume mutation
and elitism condition.

## Example initial population
- The initial population consists of binary input vectors representing
one quadrant of the actual topology, excluding the top and bottom
boundary elements.
- An example of an initial population with a volume fraction of 40%
is provided in `/matlab/GA_IP.txt`.
- The 'GA_IP.txt' file must be in the same directory as the Matlab
routines.

## Finite element software
- **This topology optimization routine uses FEAP as the finite element
software**, http://projects.ce.berkeley.edu/feap/
- The FEAP executable must be in the same directory as the Matlab
routines.

## Steps
Assuming your working directory is `/matlab`, the steps for using this
routine are:
1. Open the File `TopOptGA.m`.
2. Input the required information including:
    1. General input, `line 27` to `line 47`:
        - Geometric dimensions of the domain, including boundary elements.
		- Number of finite elements used to discretize the domain.
		- Material properties.
		- Reference load and maximum displacement.
    2. Genetic algorithm input, `line 51 to `line 77`:
		- Number of topologies.
		- Number of generations.
		- Volume fraction.
		- Connectivity condition.
		- Repair procedure.
	3. Genetic algorithm operators input, `line 79` to `line 90`:
		- Number of elite topologies, must be an even number.
		- Number of elements selected for crossover, 50% of total.
		- Probability of mutation and mutation area.
	4. Plot parameters input, `line 92` to `line 96`:
		- Selected generations for post-processing.
		- Plot configuration.
3. Load the initial population, `GA_IP.txt`.
4. Evaluate the initial population:
	1. Run repair algorithm.
	2. Construct the complete topology.
	3. Add top and bottom boundary elements.
	4. Run the image-processing-based analysis to check load path.
	5. Run shear yield FE analysis, get Vy.
	6. Run shear buckling FE analysis, get Vb.
	7. Compute and save the objective function value.
	8. Save the current topology.
5. Perform the genetic algorithm:
	1. Run the selection operation.
	2. Run the crossover operation.
	3. Run the mutation operation.
	4. Run the elitism condition.
	5. Evaluate the current population:
		- Run repair algorithm.
		- Construct the complete topology.
		- Add top and bottom boundary elements.
		- Run the image-processing-based analysis to check loadpath.
		- Run shear yield FE analysis, get Vy.
		- Run shear buckling FE analysis, get Vb.
		- Compute and save the objective function value.
		- Save the current topology.
	6. Check stop condition.
6. Post-processing:
	1. Plot Vb vs Vy and other plots using the generated information.
	2. Plot the fraction of admissible topologies.
7. Save the files as `.mat files`:
	1. Save the date as `GA_Data_400.mat`:
		- This file contains all the numerical information for each
        topology at each generation.
	2. Save the generated topologies as 'GA_Shapes_400.mat':
		- This file contains all the topologies at each generation.
		- The topologies are binary vectors representing one quadrant
        of the actual topology.


