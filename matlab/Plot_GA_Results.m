%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Controlling Out-of-Plane Buckling in Shear-Acting Structural Fuses
%%%%%% Through Topology Optimization
%%%%%% Javier A. Avecillas; Matthew R. Eatherton
%%%%%% Department of Civil and Environmental Engineering, Virginia Tech
%%%%%% Version 1.0 - Last update: 07/09/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% PLOT GENETIC ALGORITHM RESULTS %%%%%%%%%%%%%%%%%%%%%%
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
% 'shapes'        3D array with all topologies at each generation
% 'nx' and 'ny'   Number of elements in the x and y direction - Full domain
% 'nx_SD'         Number of elements in the x - Quadrant domain
% 'ny_SD'         Number of elements in the y - Quadrant domain
% 'gnts_vect'     Selected generations to be analized
% 'plot_conf'     Subplot configuration - #plots_x by #plots_y

function [ f1, f2, f3, f4 ] = Plot_GA_Results( data, shapes, nx_SD, ny_SD, nx, ny, gnts_vect, plot_conf )
%% Plost V_b vs V_y for the selected generations
f1 = figure('Name','V_buckle vs V_yield');
for i_plot = 1:length(gnts_vect)
    subplot(plot_conf(1),plot_conf(2),i_plot);
    fplot(@(x) 1.00*(x),'--k','Linewidth',1); hold on
    fplot(@(x) 0.50*(x),'--m','Linewidth',1); hold on
    fplot(@(x) 0.25*(x),'--r','Linewidth',1); hold on
    for j = 1:4
        viscircles([0,0],(j/4)*max(max(max(data(1:2,:,:)))),'Color','k','LineStyle',':','LineWidth',1);
        hold on
    end
    scatter(data(2,:,gnts_vect(i_plot)),data(1,:,gnts_vect(i_plot)),5,'b','filled');
    hold on
    
    [min_OF, idx_opt] = min(data(10,:,gnts_vect(i_plot)));
    scatter(data(2,idx_opt,gnts_vect(i_plot)),data(1,idx_opt,gnts_vect(i_plot)),20,'r','filled')

    xlim([0 1.1*max(max(data(2,:,:)))])
    ylim([0 1.1*max(max(data(1,:,:)))])
    legend('V_{Y}=1.00V_{B}','V_{Y}=0.50V_{B}','V_{Y}=0.25V_{B}','Location','Northwest','Orientation','vertical')

    set(gca,'FontName', 'Times New Roman');
    title(['Point Cloud - Generation ' num2str(gnts_vect(i_plot))])
    xlabel('Shear Buckling Load, V_B (kips)','FontName', 'Times New Roman')
    ylabel('Shear Yielding Load, V_Y (kips)','FontName', 'Times New Roman')
end

%% Plots V_b vs V_y vs OF for the selected generations
f2 = figure('Name','Objective Function');
for i_plot = 1:length(gnts_vect)
    subplot(plot_conf(1),plot_conf(2),i_plot);
    scatter3(data(2,:,gnts_vect(i_plot)),data(1,:,gnts_vect(i_plot)),data(10,:,gnts_vect(i_plot)),5,'b','filled'); hold on
    
    set(gca,'FontName', 'Times New Roman');
    title(['Objective Function - Generation ' num2str(gnts_vect(i_plot))])
    xlabel('V_B (kips)','FontName', 'Times New Roman')
    ylabel('V_Y (kips)','FontName', 'Times New Roman')
    zlabel('OF','FontName', 'Times New Roman')
    xlim([0 1.1*max(max(data(2,:,:)))])
    ylim([0 1.1*max(max(data(1,:,:)))])
    zlim([0 1.0])
    
    [min_OF, idx_opt] = min(data(10,:,gnts_vect(i_plot)));
    scatter3(data(2,idx_opt,gnts_vect(i_plot)),data(1,idx_opt,gnts_vect(i_plot)),data(10,idx_opt,gnts_vect(i_plot)),20,'r','filled')
end

%% Plot best topology at each selected generation
f3 = figure('Name','Best Topologies');
for i_plot = 1:length(gnts_vect)
    subplot(plot_conf(1),plot_conf(2),i_plot);
    [min_OF, idx_opt] = min(data(10,:,gnts_vect(i_plot)));
    opt_top = shapes(:,idx_opt,gnts_vect(i_plot));
    % Apply simmetry
    symm = 3;
    [ S_opt_top ] = Symmetry_Input( opt_top, nx_SD, ny_SD, symm );
    % Apply boundary elements
    B_S_opt_top = vertcat(ones(2*nx_SD,1),S_opt_top,ones(2*nx_SD,1));
    % From binary vector to binary matrix
    B_S_opt_top = (reshape(B_S_opt_top,[nx,ny]))';
    B_S_opt_top = flip(B_S_opt_top,1);
    BW = 1-B_S_opt_top;
    S_BW = imresize(BW,1);
    imshow(S_BW)

    set(gca,'FontName', 'Times New Roman');
    title(['Optimum Topology - Generation ' num2str(gnts_vect(i_plot))])
    xlabel(['Objective Function = ' num2str(min_OF)])
end

%% Plot OF values through generations
f4 = figure('Name','OF Evoution');
n_gnts = size(data,3);
OF_History = zeros(1,n_gnts);
for i = 1:n_gnts
    [min_OF, idx_opt] = min(data(10,:,i));
    OF_History(1,i) = min_OF;
end
plot(linspace(1,n_gnts,n_gnts),OF_History)
xlim([1 n_gnts])
set(gca,'FontName', 'Times New Roman');
title('OF Evolution')
xlabel('Generation','FontName', 'Times New Roman')
ylabel('Objective Funciton','FontName', 'Times New Roman')

end