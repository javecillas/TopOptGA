%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Controlling Out-of-Plane Buckling in Shear-Acting Structural Fuses
%%%%%% Through Topology Optimization
%%%%%% Javier A. Avecillas; Matthew R. Eatherton
%%%%%% Department of Civil and Environmental Engineering, Virginia Tech
%%%%%% Version 1.0 - Last update: 07/09/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% FRACTION OF ADMISSIBLE TOPOLOGIES %%%%%%%%%%%%%%%%%%%%
% 'data'          Information about each topology at each generation
%                 data(1,i) -> Vy
%                 data(2,i) -> Vb
%                 data(3,i) -> sqrt(Vy^2+Vb^2)
%                 data(4,i) -> Number of active components
%                 data(5,i) -> Load path condition 1-yes 0-no
%                 data(6,i) -> Number of disconnected components
%                 data(7,i) -> Area of disconnected components
%                 data(8,i) -> Min. area of disconnected components
%                 data(9,i) -> Min. area of internal holes
%                 data(10,i) -> Objective function value

function [ f1 ] = Plot_Feasible( data )
%% Step 1
% Get the total number of generations and individuals per generation
gnts = size(data,3);
nind = size(data,2);

%% Step 2
% Get the fraction of the total population corresponding to admissible
% topologies
feas = zeros(1,gnts);
for c_gnt = 1:gnts
    feas(1,c_gnt) = (1/nind)*(sum(data(5,:,c_gnt) == 1 & data(6,:,c_gnt) == 0));
end

%% Step 3
% Plot
f1 = figure('Name','Feasible Individuals Evolution');
plot(linspace(1,gnts,gnts),feas);
xlim([1 gnts])
ylim([0 1])

set(gca,'FontName', 'Times New Roman');
title('Feasible Individuals Evolution')
xlabel('Generation','FontName', 'Times New Roman')
ylabel('Admissible topologies','FontName', 'Times New Roman')

end