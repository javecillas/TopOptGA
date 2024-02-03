%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Controlling Out-of-Plane Buckling in Shear-Acting Structural Fuses
%%%%%% Through Topology Optimization
%%%%%% Javier A. Avecillas; Matthew R. Eatherton
%%%%%% Department of Civil and Environmental Engineering, Virginia Tech
%%%%%% Version 1.0 - Last update: 07/09/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% SHEAR YIELDING ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% This Matlab code needs FEAP (Taylor 2013) to calculate Vy
%%%%%% The FE analysis only uses one quadrant of the topology
% 'input'         2D binary array representing the repaired topology
%                 including the top and bottom boundary elements
% 'width'         Width of the domain
% 'height'        Height of the domain
% 'thickness'     Thickness of the domain
% 'nx' and 'ny'   Number of elements in the x and y direction
% 'E'             Material Young's modulus
% 'Fy'            Material yield strength
% 'FEAP_ref_disp' Maximum applied displacement for shear yield analysis 
% 'fetopo'        Elemental connectivity matrix                          %
% 'fecoord'       Nodal global coordinates 

function [ FEAP_shear_yielding ] = FEAP_Shear_Yielding( Input, width, height, thickness, nx, ny, E, Fy, FEAP_ref_disp, fetopo, fecoord )
%% Matrix array with node numbers
% Counter clockwise node numbering
nod_numb = reshape(linspace(1,(nx+1)*(ny+1),(nx+1)*(ny+1)),nx+1,ny+1)';

%% Recall nodes coordinates
% Counter clockwise mesh nodes coordinates
aux_node = zeros(4,(nx+1)*(ny+1));
for m = 1:nod_numb(end,end)
    aux_node(1,m) = m;
    aux_node(3,m) = fecoord(m,1);
    aux_node(4,m) = fecoord(m,2);
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

%% Essential boundary conditions - Doubly simmetry
% Bottom left corner & Top and Bottom boundary for all nodes
% Left and Right boundary for all nodes besides Top and Bottom corners
% For symmetrical FEA - Bottom Left Corner Essential Boundary Condition
sym_boun_cor = 1;
sym_boun_cor(find(ismember(sym_boun_cor,aux_ele_edited(4:end,:)) == 0)) = [];
if isempty(sym_boun_cor) == 1
else
    sym_boun_cor = vertcat(sym_boun_cor,zeros(1,size(sym_boun_cor,2)));
    sym_boun_cor = vertcat(sym_boun_cor,repmat([1 1]',1,size(sym_boun_cor,2)));
end
% For symmetrical FEA - Bottom Edge Essential Boundary Condition
sym_boun_bot = linspace(2,(nx+1),(nx));
sym_boun_bot(find(ismember(sym_boun_bot,aux_ele_edited(4:end,:)) == 0)) = [];
if isempty(sym_boun_bot) == 1
else
    sym_boun_bot = vertcat(sym_boun_bot,zeros(1,size(sym_boun_bot,2)));
    sym_boun_bot = vertcat(sym_boun_bot,repmat([1 0]',1,size(sym_boun_bot,2)));
end
% For symmetrical FEA - Right Edge Essential Boundary Condition
sym_boun_rig = linspace(2*(nx+1),(nx+1)+(ny-1)*(nx+1),(ny-1));
sym_boun_rig(find(ismember(sym_boun_rig,aux_ele_edited(4:end,:)) == 0)) = [];
if isempty(sym_boun_rig) == 1
else
    sym_boun_rig = vertcat(sym_boun_rig,zeros(1,size(sym_boun_rig,2)));
    sym_boun_rig = vertcat(sym_boun_rig,repmat([0 0]',1,size(sym_boun_rig,2)));
end
% For symmetrical FEA - Top Edge Essential Boundary Condition
sym_boun_top = linspace(((nx+1)*(ny+1)-nx),((nx+1)*(ny+1)),(nx+1));
sym_boun_top(find(ismember(sym_boun_top,aux_ele_edited(4:end,:)) == 0)) = [];
if isempty(sym_boun_top) == 1
else
    sym_boun_top = vertcat(sym_boun_top,zeros(1,size(sym_boun_top,2)));
    sym_boun_top = vertcat(sym_boun_top,repmat([0 1]',1,size(sym_boun_top,2)));
end
% For symmetrical FEA - Left Edge Essential Boundary Condition
sym_boun_lef = linspace(nx+2,(1+(ny-1)*(nx+1)),ny-1);
sym_boun_lef(find(ismember(sym_boun_lef,aux_ele_edited(4:end,:)) == 0)) = [];
if isempty(sym_boun_lef) == 1
else
    sym_boun_lef = vertcat(sym_boun_lef,zeros(1,size(sym_boun_lef,2)));
    sym_boun_lef = vertcat(sym_boun_lef,repmat([0 1]',1,size(sym_boun_lef,2)));
end
% For symmetrical FEA - Essencial Boundary Conditions - Summary
aux_fix = [sym_boun_cor sym_boun_bot sym_boun_top sym_boun_lef];

%% Fix unconnected nodes
% Look for nodes which are not part of the current mesh and current bounds
for i = 1:(nx+1)*(ny+1)
    if ismember(i,aux_ele_edited(4:7,:)) == 0 && ismember(i,aux_fix(1,:)) == 0
        aux_fix = [aux_fix [i 0 1 1]'];
    else
    end
end

%% Asign loaded nodes
% The load is an unit load in the x-dir distributed along the top edge
aux_disp = zeros(4,(nx+1));
aux_disp(1,:) = linspace((nx+1)*(ny+1)-nx,(nx+1)*(ny+1),(nx+1));
aux_disp(2,:) = 0;
% Set Reference Load for shear yield analysis
% distributed among all the Top edge's nodes
aux_disp(3,:) = 1;

%% Print FEAP input txt file
fid = fopen('IFEAP_shear_yielding.txt','w');
fprintf(fid,'FEAP * * 3-D Shell Problem\n');
formatSpec =('%5.0f %5.0f %5.0f %1.0f %1.0f %1.0f\n');
fprintf(fid,formatSpec, [0; 0; 0; 2; 2; 4]);
fprintf(fid,'\n');
%%%%%%%%%%%%%%%%%%%%%%%% ELEMENT PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'NOPRint\n');
fprintf(fid,'MATErial,1\n');
fprintf(fid,'SOLId\n');
fprintf(fid,'PLANe STREss\n');

FEAP_material_properties = ['ELAStic ISOTropic' ' ' num2str(E) ' ' '0.3\n'];
fprintf(fid,FEAP_material_properties);

FEAP_mises_properties = ['PLAStic MISEs' ' ' num2str(Fy) ' ' num2str(Fy) ' ' num2str(0.0) ' ' '\n'];
fprintf(fid,FEAP_mises_properties);

FEAP_hardening_properties = ['PLAStic HARDening' ' ' num2str(0.0) ' ' num2str(0.01*E) ' ' '\n'];
fprintf(fid,FEAP_hardening_properties);

fprintf(fid,'QUADrature order 2 2\n');

fprintf(fid,'\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nodes for elements
fprintf(fid,'NOPRint\n');
fprintf(fid,'COORdinates\n');
formatSpec =('%5.0f %1.0f %5.5f %5.5f\n');
fprintf(fid,formatSpec, aux_node);
fprintf(fid,'\n');

% 2D Plane-Stress elements
fprintf(fid,'ELEMents\n');
formatSpec =('%5.0f %1.0f %5.0f %5.0f %5.0f %5.0f %5.0f\n');
fprintf(fid,formatSpec, aux_ele_edited);
fprintf(fid,'\n');

% Equal DOF at Top Edge's nodes
fprintf(fid,'LINK\n');
FEAP_link = [num2str((nx+1)*(ny+1)) ',' num2str((nx+1)*(ny+1)-nx) ',' '0' ',' '1' ',' '1' ',' '0' '\n'];
fprintf(fid,FEAP_link);
FEAP_link = [num2str((nx+1)*(ny+1)) ',' num2str((nx+1)*(ny+1)-nx) ',' ',' ',' '1' ',' '0' '\n'];
fprintf(fid,FEAP_link);
fprintf(fid,'\n');

% Boundary conditions
fprintf(fid,'BOUNdary restraints\n');
formatSpec =('%5.0f %1.0f %1.0f %1.0f\n');
fprintf(fid,formatSpec, aux_fix);
fprintf(fid,'\n');

% Applied reference load
fprintf(fid,'FORCes\n');
formatSpec =('%5.0f %1.0f %5.5f %5.5f\n');
fprintf(fid,formatSpec, aux_disp);
fprintf(fid,'\n');

fprintf(fid,'END\n');
fprintf(fid,'\n');

fprintf(fid,'BATCh\n');

fprintf(fid,'ARCLength,,5\n');
fprintf(fid,'DT,,0.1\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'LOOP,,100\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'TIME\n');
fprintf(fid,'LOOP,,10\n');
fprintf(fid,'TANG,,1\n');
fprintf(fid,'NEXT\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
formatSpec = 'DISP,,%5.0f\n';
fprintf(fid,formatSpec,(nx+1)*(ny+1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'REAC,coor,2,0.0\n');
fprintf(fid,'NEXT\n');
fprintf(fid,'END\n');

fprintf(fid,'0\n');

formatSpec = '%5.0f %1.0f %1.5f\n';
fprintf(fid,formatSpec,(nx+1)*(ny+1),1,FEAP_ref_disp/100);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf(fid,'INTER\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'STOP\n');
fprintf(fid,'\n');

fclose(fid);

%% Run FEAP84
! feap84.exe -iIFEAP_shear_yielding.txt -oOFEAP_shear_yielding.txt

%% Post Processing
% Calculate the shear yield load
FEAP_Out_File = fileread('OFEAP_shear_yielding.txt');
FEAP_Out_File = regexp(FEAP_Out_File,'\n','split')';

txt_index_c = strfind(FEAP_Out_File, 'N o d a l   D i s p l a c e m e n t s');
txt_index_disp = find(not(cellfun('isempty', txt_index_c)));

txt_index_c = strfind(FEAP_Out_File, 'N o d a l    R e a c t i o n s');
txt_index_reac = find(not(cellfun('isempty', txt_index_c)));

FEAP_shear_yielding_disp  = str2num(cell2mat(FEAP_Out_File(txt_index_disp+4,1)));

FEAP_shear_yielding_forc = [];
for i = 3:3+(size([sym_boun_cor sym_boun_bot],2)-1)
    FEAP_shear_yielding_forc  = [FEAP_shear_yielding_forc str2num(cell2mat(FEAP_Out_File(txt_index_reac+i,1)))];
end

FEAP_shear_yielding_forc_plot = zeros(length(FEAP_shear_yielding_disp),1);
for i = 1:length(FEAP_shear_yielding_disp)
    sum_fx = 0;
    for j = 0:(size([sym_boun_cor sym_boun_bot],2)-1)
        sum_fx = sum_fx+FEAP_shear_yielding_forc(i,2+3*j);
    end
    FEAP_shear_yielding_forc_plot(i) = -1*2*thickness*sum_fx;
end

% Compute the shear yield load
FEAP_shear_yielding = 0;
FEAP_ini_stiffness = (FEAP_shear_yielding_forc_plot(2,1)-FEAP_shear_yielding_forc_plot(1,1))/(FEAP_shear_yielding_disp(2,4)-FEAP_shear_yielding_disp(1,4));
if FEAP_ini_stiffness == 0
    FEAP_shear_yielding = 0;
    return
else
    for i = 2:length(FEAP_shear_yielding_disp)-1
        FEAP_current_stiffness = (FEAP_shear_yielding_forc_plot(i+1,1)-FEAP_shear_yielding_forc_plot(i,1))/(FEAP_shear_yielding_disp(i+1,4)-FEAP_shear_yielding_disp(i,4));
        if FEAP_current_stiffness < 0.25*FEAP_ini_stiffness
            FEAP_shear_yielding = FEAP_shear_yielding_forc_plot(i+1);
            break
        else
        end
    end
end

end

