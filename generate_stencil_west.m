
% This code generates a 2-D FVM stencil for inner nodes.
% The faces of the cells do not have to be aligned to a 
% cartesian grid on any side, i.e. a cell can be any convex quadrilateral

% The code is devided in two parts. The first part generates the stencil
% and prints it to a file given in variable 'target file' after the 
% expression 'Start_stecil' (case sensivtive). The second part
% replaces the scalar expressions with matrix expression and allows to
% calculate the stencil for all inner nodes at once. This accelerates the
% calculation of the system Matrix a lot. 

% Either way all distances and surfaces (areas) have to be provided by you.

% Written by Thomas Runte


clear all 
clc

% Print results to file 
target_file = 'build_west.m';

fclose(fopen(target_file, 'w'))

%% First part
% Generate stencil with variables names introduced in Camilo F. Silva's
% Course "Numerical Thermo Fluid Dynamics"

% Initialize symbolic variables for distances,
% areas, lambdas, Temperatures



% Around sepsilon
syms dy_S_Se dy_Se_e dy_e_P dy_P_S  real
syms dx_S_Se dx_Se_e dx_e_P dx_P_S  real
            
% Around nepsilon
syms dy_P_e dy_e_Ne dy_Ne_N dy_N_P  real            
syms dx_P_e dx_e_Ne dx_Ne_N dx_N_P  real
            
% Around e
syms dy_s_sE dy_sE_nE dy_nE_n dy_n_s real
syms dx_s_sE dx_sE_nE dx_nE_n dx_n_s real
            
% Around P
syms dy_s_se dy_se_ne dy_ne_n dy_n_s  real            
syms dx_s_se dx_se_ne dx_ne_n dx_n_s  real


% Areas
syms S_e S_sepsilon S_nepsilon  S_epsilon  real

% Temperatures
syms  T_E T_W T_S T_N T_NW T_NE T_SW T_SE T_P real

% inner Temperatures
syms T_sw T_se T_ne T_nw real

% Temperatures at the boundaries
syms T_Nepsilon T_Sepsilon T_epsilon real

% Define inner Temperatures as interpolation of outer Temperatures
T_sw=(T_SW+T_S+T_P+T_W)/4;
T_se=(T_S+T_SE+T_E+T_P)/4;
T_ne=(T_P+T_E+T_NE+T_N)/4;
T_nw=(T_W+T_P+T_N+T_NW)/4;

T_s    = 0.5*T_S + 0.5*T_P;
T_n    = 0.5*T_P + 0.5*T_N;
T_Nepsilon = 0.75*T_N + 0.25*T_NE;
T_Sepsilon =  0.75*T_S + 0.25*T_SE;
T_epsilon = 0.75*T_P + 0.25*T_E;

% Gradients (Greens theorem)

dTdx_sepsilon=  (dy_S_Se*T_Sepsilon + dy_Se_e*T_se + dy_e_P*T_epsilon + dy_P_S*T_s)/S_sepsilon;
dTdy_sepsilon= -(dx_S_Se*T_Sepsilon + dx_Se_e*T_se + dx_e_P*T_epsilon + dx_P_S*T_s)/S_sepsilon;


dTdx_nepsilon=  (dy_P_e*T_epsilon + dy_e_Ne*T_ne + dy_Ne_N*T_Nepsilon + dy_N_P*T_n)/S_nepsilon;
dTdy_nepsilon= -(dx_P_e*T_epsilon + dx_e_Ne*T_ne + dx_Ne_N*T_Nepsilon + dx_N_P*T_n)/S_nepsilon;


dTdx_e=   (dy_s_sE*T_se + dy_sE_nE*T_E + dy_nE_n*T_ne + dy_n_s*T_P)/S_e;
dTdy_e=  -(dx_s_sE*T_se + dx_sE_nE*T_E + dx_nE_n*T_ne + dx_n_s*T_P)/S_e;


% Build hole stecil acounting for quadratic lambda like in Helmholtz 

  DDT =( dy_s_se*dTdx_sepsilon - dx_s_se*dTdy_sepsilon ...
       + dy_se_ne*dTdx_e - dx_se_ne*dTdy_e ...
       + dy_ne_n*dTdx_nepsilon     - dx_ne_n*dTdy_nepsilon)/S_epsilon;

 
% Make Temperatur vector 
T=[T_E; T_W; T_S; T_N; T_NW; T_NE; T_SW; T_SE; T_P];


% NOTE: The Jacobian function is a "smart" way of factorizing the
% expresssion DDT with respect the temperature vector

stecil=jacobian(DDT,T);

% Find position in file
fileID2 = fopen(target_file, 'r+');

        fprintf(fileID2,'\n\n');
        fprintf(fileID2,'%% Nomenclature:\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%%     N(i-1,j) -  Ne     NE(i-1,j+1)\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%%       |         |         |\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%%      n ------  ne - - - nE\n');
        fprintf(fileID2,'%%                 |                 |\n');
        fprintf(fileID2,'%%       |         |        |        |       |\n');
        fprintf(fileID2,'%%                 |                 |\n');
        fprintf(fileID2,'%%   P (i,j) - -   e - -  E (i,j+1)\n');
        fprintf(fileID2,'%%                 |                 |\n');
        fprintf(fileID2,'%%       |         |        |        |       |\n');
        fprintf(fileID2,'%%                 |                 |\n');
        fprintf(fileID2,'%%       s ------ se - - - sE\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%%       |         |        |        |\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%%     S(i+1,j)  - Se      SE(i+1,j+1)\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%% Indexing of stecil: \n\n');
        fprintf(fileID2,'%%   D_1 - D2 \n');
        fprintf(fileID2,'%%     |     | \n\');
        fprintf(fileID2,'%%   D_0 - D3\n');
        fprintf(fileID2,'%%     |     | \n\');
        fprintf(fileID2,'%%    D1 - D4\n\n');

        
        %fprintf(fileID2,'lambda=boundary.lambda; \n\n');
      
        fprintf(fileID2,'%% Stecil \n\n');
        fprintf(fileID2,'%% East \n');
        fprintf(fileID2,'D3=%s; \n\n',char(stecil(1)));
   
        fprintf(fileID2,'%% West \n');
        fprintf(fileID2,'D_3=%s; \n\n',char(stecil(2)));

        fprintf(fileID2,'%% South \n');
        fprintf(fileID2,'D1=%s; \n\n',char(stecil(3)));

        fprintf(fileID2,'%% North \n');
        fprintf(fileID2,'D_1=%s; \n\n',char(stecil(4)));

        fprintf(fileID2,'%% NW \n');
        fprintf(fileID2,'D_4=%s; \n\n',char(stecil(5)));

        fprintf(fileID2,'%% NE \n');
        fprintf(fileID2,'D2=%s; \n\n',char(stecil(6)));

        fprintf(fileID2,'%% SW \n');
        fprintf(fileID2,'D_2=%s; \n\n',char(stecil(7)));

        fprintf(fileID2,'%% SE \n');
        fprintf(fileID2,'D4=%s; \n\n',char(stecil(8)));

        fprintf(fileID2,'%% P \n');
        fprintf(fileID2,'D0=%s; \n\n',char(stecil(9)));

fclose(fileID2);
