%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Controlling Out-of-Plane Buckling in Shear-Acting Structural Fuses
%%%%%% Through Topology Optimization
%%%%%% Javier A. Avecillas; Matthew R. Eatherton
%%%%%% Department of Civil and Environmental Engineering, Virginia Tech
%%%%%% Version 1.0 - Last update: 07/09/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% LINEAR BUCKLING ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% This Matlab code needs FEAP (Taylor 2013) to calculate Vb
%%%%%% The FE analysis uses the full geometry of the topology
% 'input'         2D binary array representing the repaired topology
%                 including the top and bottom boundary elements
% 'thickness'     Thickness of the domain
% 'nx' and 'ny'   Number of elements in the x and y direction
% 'E'             Material Young's modulus
% 'FEAP_ref_load' Reference load for eigen-buckling analysis 
% 'fetopo'        Elemental connectivity matrix                          %
% 'fecoord'       Nodal global coordinates 

function [ FEAP_shear_buckling ] = FEAP_Shear_Buckling( Input, thickness, nx, ny, E, FEAP_ref_load, fetopo, fecoord )
%% Matrix array with node numbers
% Counter Clockwise Node Numbering
nod_numb = reshape(linspace(1,(nx+1)*(ny+1),(nx+1)*(ny+1)),nx+1,ny+1)';

%% Recall nodes coordinates
% Counter clockwise mesh nodes coordinates
% Initial out-of-plane imperfection
z_0 = 0;
aux_node = zeros(5,(nx+1)*(ny+1));
for m = 1:nod_numb(end,end)
    aux_node(1,m) = m;
    aux_node(3,m) = fecoord(m,1);
    aux_node(4,m) = fecoord(m,2);
    aux_node(5,m) = z_0;
end

%% Recall element connectivity
aux_ele = zeros(7,nx*ny);
aux_ele(1,:) = linspace(1,(nx*ny),(nx*ny));
% Nodal increment
aux_ele(2,:) = 0;
aux_ele(3,:) = 1;
aux_ele(4:end,:) = fetopo';
% Duplicated mesh - Used to delete void elements and renumbering
aux_ele_edited = aux_ele;

%% Delete void elements
% Based on the duplicated mesh
erase = find(Input(:,1) == 0);
aux_ele_edited(:,erase) = [];
% Renumbering the material elements
aux_ele_edited(1,:) = linspace(1,size(aux_ele_edited,2),size(aux_ele_edited,2));

%% Essential boundary conditions
% Bottom edge essential boundary condition
boun_bot = zeros(8,(nx+1));
boun_bot(1,:) = linspace(1,(nx+1),(nx+1));
boun_bot(3:8,:) = repmat([1 1 1 1 1 1]',1,size(boun_bot,2));
% Rigth edge essential boundary condition
boun_rig = zeros(8,(ny+1));
boun_rig(1,:) = linspace((nx+1),(nx+1)+ny*(nx+1),(ny+1));
boun_rig(3:8,:) = repmat([0 0 1 0 0 0]',1,size(boun_rig,2));
% Top edge essential boundary condition
boun_top = zeros(8,(nx+1));
boun_top(1,:) = linspace(((nx+1)*(ny+1)-nx),((nx+1)*(ny+1)),(nx+1));
boun_top(3:8,:) = repmat([0 1 1 1 1 1]',1,size(boun_top,2));
% Left edge essential boundary condition
boun_lef = zeros(8,(ny+1));
boun_lef = linspace(1,(1+ny*(nx+1)),(ny+1));
boun_lef(3:8,:) = repmat([0 0 1 0 0 0]',1,size(boun_lef,2));
% Essencial Boundary Conditions - Summary
aux_fix = [boun_bot boun_top];

%% Fix unconnected nodes
% Look for nodes which are not part of the current mesh 
for i = 1:(nx+1)*(ny+1)
    if ismember(i,aux_ele_edited(4:7,:)) == 0
        aux_fix = [aux_fix [i 0 1 1 1 1 1 1]'];
    else
    end
end

%% Asign loaded nodes
% The load is distributed among the top edge nodes
aux_force = zeros(8,(nx+1));
aux_force(1,:) = linspace((nx+1)*(ny+1)-nx,(nx+1)*(ny+1),(nx+1));
aux_force(2,:) = 0;
% Set Reference Load for linear buckling analysis
% distributed among all the Top edge's nodes
aux_force(3,1) = FEAP_ref_load/(2*nx);
aux_force(3,2:end-1) = FEAP_ref_load/nx;
aux_force(3,end) = FEAP_ref_load/(2*nx);

%% Print FEAP input txt file
fid = fopen('IFEAP_shear_buckling.txt','w');
fprintf(fid,'FEAP * * 3-D Shell Problem\n');
formatSpec =('%5.0f %5.0f %5.0f %1.0f %1.0f %1.0f\n');
fprintf(fid,formatSpec, [0; 0; 0; 3; 6; 4]);
fprintf(fid,'\n');
%%%%%%%%%%%%%%%%%%%%%%%%% SHELL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'NOPRint\n');
fprintf(fid,'MATErial,1\n');
fprintf(fid,'SHELl\n');

FEAP_material_properties = ['ELAStic ISOTropic' ' ' num2str(E) ' ' '0.3\n'];
fprintf(fid,FEAP_material_properties);

FEAP_thickness = ['THICk,,' ' ' num2str(thickness) ' ' '5/6 5\n'];
fprintf(fid,FEAP_thickness);

fprintf(fid,'\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nodes for shell elements
fprintf(fid,'NOPRint\n');
fprintf(fid,'COORdinates\n');
formatSpec =('%5.0f %1.0f %5.5f %5.5f %5.5f\n');
fprintf(fid,formatSpec, aux_node);
fprintf(fid,'\n');

% Shell elements
fprintf(fid,'ELEMents\n');
formatSpec =('%5.0f %1.0f %5.0f %5.0f %5.0f %5.0f %5.0f\n');
fprintf(fid,formatSpec, aux_ele_edited);
fprintf(fid,'\n');

% Boundary conditions
fprintf(fid,'BOUNdary restraints\n');
formatSpec =('%5.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f\n');
fprintf(fid,formatSpec, aux_fix);
fprintf(fid,'\n');

% Applied reference load
fprintf(fid,'FORCe\n');
formatSpec =('%5.0f %1.0f %5.5f %5.5f %5.5f %5.5f %5.5f %5.5f\n');
fprintf(fid,formatSpec, aux_force);
fprintf(fid,'\n');

fprintf(fid,'END\n');
fprintf(fid,'\n');

fprintf(fid,'BATCh\n');

fprintf(fid,'LOOP,,100\n');
fprintf(fid,'TANGent,,1\n');
fprintf(fid,'NEXT\n');
%fprintf(fid,'OUTP,TANGent\n');

fprintf(fid,'GEOMetric\n');
%%%%%%%%%%%%%%%%%%%%%%%% SUBSPACE ALGORITHM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'SUBS,,2,20,1e-12,50\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% PLOT EIGEN BUCKLING MODES %%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf(fid,'EIGV,all,1\n');
% fprintf(fid,'EIGV,all,2\n');
% fprintf(fid,'EIGV,all,3\n');
% fprintf(fid,'EIGV,all,4\n');
% 
% fprintf(fid,'PLOT Frame 1\n');
% fprintf(fid,'PLOT EIGV,1,-1,3\n');
% fprintf(fid,'PLOT Frame 2\n');
% fprintf(fid,'PLOT EIGV,2,-1,3\n');
% fprintf(fid,'PLOT Frame 3\n');
% fprintf(fid,'PLOT EIGV,3,-1,3\n');
% fprintf(fid,'PLOT Frame 4\n');
% fprintf(fid,'PLOT EIGV,4,-1,3\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'END\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf(fid,'INTEractive\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'STOP\n');
fprintf(fid,'\n');

fclose(fid);

%% Run FEAP84
! feap84.exe -iIFEAP_shear_buckling.txt -oOFEAP_shear_buckling.txt

%% Post Processing
% Calculate the shear buckling load
FEAP_Out_File = fileread('OFEAP_shear_buckling.txt');
FEAP_Out_File = regexp(FEAP_Out_File,'\n','split')';

txt_index_c = strfind(FEAP_Out_File, 'SUBSPACE: Current eigenvalues');
txt_index_subspace = find(not(cellfun('isempty', txt_index_c)));

txt_index_c = strfind(FEAP_Out_File, 'SUBSPACE: Current residuals');
txt_index_residuals = find(not(cellfun('isempty', txt_index_c)));

txt_subspace_ite = (regexp(cell2mat(FEAP_Out_File(txt_index_subspace,1)),'\d*','Match'));
subspace_ite = str2num(txt_subspace_ite{1,1});

if subspace_ite == 50;
    FEAP_shear_buckling  = 0;
    return
else
    FEAP_lambda_buckling = [];
    for c_lambda = txt_index_subspace+1:txt_index_residuals-2
        FEAP_lambda_buckling = [FEAP_lambda_buckling; (str2num(cell2mat(FEAP_Out_File(c_lambda,1))))'];
    end

    FEAP_residuals = [];
    for c_res = txt_index_residuals+1:(txt_index_residuals+1)+((txt_index_residuals-2)-(txt_index_subspace+1))
        FEAP_residuals = [FEAP_residuals; (str2num(cell2mat(FEAP_Out_File(c_res,1))))'];
    end

    FEAP_shear_buckling  = 0;
    for i = 1:numel(FEAP_lambda_buckling)
    if abs(FEAP_residuals(i)) < 1.0*10^-10
        FEAP_shear_buckling  = abs(FEAP_ref_load*FEAP_lambda_buckling(i));
        break
    else
    end
    end
end

end

