
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
target_file = 'build_east.m';

fclose(fopen(target_file, 'w'))

%% First part
% Generate stencil with variables names introduced in Camilo F. Silva's
% Course "Numerical Thermo Fluid Dynamics"

% Initialize symbolic variables for distances,
% areas, lambdas, Temperatures



% Around somega
syms dy_Sw_S dy_S_P dy_P_w dy_w_Sw real
syms dx_Sw_S dx_S_P dx_P_w dx_w_Sw real
            
% Around nomega
syms dy_w_P dy_P_N dy_N_Nw dy_Nw_w  real            
syms dx_w_P dx_P_N dx_N_Nw dx_Nw_w  real
            
% Around w
syms dy_sW_s dy_s_n dy_n_nW dy_nW_sW real
syms dx_sW_s dx_s_n dx_n_nW dx_nW_sW real
            
% Around P
syms dy_sw_s dy_s_n dy_n_nw dy_nw_sw  real            
syms dx_sw_s dx_s_n dx_n_nw dx_nw_sw  real


% Areas
syms S_w S_somega S_nomega  S_omega  real

% Temperatures
syms  T_E T_W T_S T_N T_NW T_NE T_SW T_SE T_P real

% inner Temperatures
syms T_sw T_se T_ne T_nw real

% Temperatures at the boundaries
syms T_Nomega T_Somega T_omega real

% Define inner Temperatures as interpolation of outer Temperatures
T_sw=(T_SW+T_S+T_P+T_W)/4;
T_se=(T_S+T_SE+T_E+T_P)/4;
T_ne=(T_P+T_E+T_NE+T_N)/4;
T_nw=(T_W+T_P+T_N+T_NW)/4;

T_s    = 0.5*T_S + 0.5*T_P;
T_n    = 0.5*T_P + 0.5*T_N;
T_Nomega = 0.75*T_N + 0.25*T_NW;
T_Somega = 0.75*T_S + 0.25*T_SW;
T_omega = 0.75*T_P + 0.25*T_W;

% Gradients (Greens theorem)

dTdx_somega=  (dy_S_P*T_s+dy_P_w*T_omega+dy_w_Sw*T_sw+dy_Sw_S*T_Somega)/S_somega;
dTdy_somega= -(dx_S_P*T_s+dx_P_w*T_omega+dx_w_Sw*T_sw+dx_Sw_S*T_Somega)/S_somega;

dTdx_nomega=  (dy_P_N*T_n+dy_N_Nw*T_Nomega+dy_Nw_w*T_nw+dy_w_P*T_omega)/S_nomega;
dTdy_nomega= -(dx_P_N*T_n+dx_N_Nw*T_Nomega+dx_Nw_w*T_nw+dx_w_P*T_omega)/S_nomega;

dTdx_w=       (dy_s_n*T_P+dy_n_nW*T_nw+dy_nW_sW*T_W+dy_sW_s*T_sw)/S_w;
dTdy_w=      -(dx_s_n*T_P+dx_n_nW*T_nw+dx_nW_sW*T_W+dx_sW_s*T_sw)/S_w;

% Build hole stecil acounting for quadratic lambda like in Helmholtz 

  DDT =(dTdx_somega*dy_sw_s-dTdy_somega*dx_sw_s ...
      + dTdx_nomega*dy_n_nw-dTdy_nomega*dx_n_nw ...
      + dTdx_w*dy_nw_sw-dTdy_w*dx_nw_sw)/S_omega ;
 
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
        fprintf(fileID2,'%%   NW=(i-1,j-1)    Nw -Nomega - N=(i-1,j)\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%%                   |            |\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%%       nW   - - - -nw- nomega - n\n');
        fprintf(fileID2,'%%                   |\n');
        fprintf(fileID2,'%%       |           |            |\n');
        fprintf(fileID2,'%%                   |\n');
        fprintf(fileID2,'%%    W=(i,j-1) - -  w - omega   P=(i,j)\n');
        fprintf(fileID2,'%%                   |\n');
        fprintf(fileID2,'%%       |           |            |\n');
        fprintf(fileID2,'%%                   |\n');
        fprintf(fileID2,'%%       sW - - -    sw- somega - s\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%%                   |            |\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%%    SW=(i+1,j-1)   Sw -Somega S=(i+1,j)\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%% Indexing of stecil: \n\n');
        fprintf(fileID2,'%%    D_4 - D_1 \n');
        fprintf(fileID2,'%%     |     |     | \n\');
        fprintf(fileID2,'%%    D_3 - D_0 \n');
        fprintf(fileID2,'%%     |     |     | \n\');
        fprintf(fileID2,'%%    D_2 -  D1 \n\n');

        
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
