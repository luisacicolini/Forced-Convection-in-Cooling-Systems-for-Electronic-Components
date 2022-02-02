
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
target_file = 'build_inner.m';

fclose(fopen(target_file, 'w'))

%% First part
% Generate stencil with variables names introduced in Camilo F. Silva's
% Course "Numerical Thermo Fluid Dynamics"

% Initialize symbolic variables for distances,
% areas, lambdas, Temperatures

% Around s
syms dy_Sw_Se dy_Se_e dy_e_w dy_w_Sw real
syms dx_Sw_Se dx_Se_e dx_e_w dx_w_Sw real

% Around e
syms dy_s_sE dy_sE_nE dy_nE_n dy_n_s real
syms dx_s_sE dx_sE_nE dx_nE_n dx_n_s real

% Around n
syms dy_w_e dy_e_Ne dy_Ne_Nw dy_Nw_w  real
syms dx_w_e dx_e_Ne dx_Ne_Nw dx_Nw_w  real

% Around w

syms dy_sW_s dy_s_n dy_n_nW dy_nW_sW real
syms dx_sW_s dx_s_n dx_n_nW dx_nW_sW real 

% Around P
syms dx_sw_se dx_se_ne dx_ne_nw dx_nw_sw real
syms dy_sw_se dy_se_ne dy_ne_nw dy_nw_sw real

% Areas
syms S_P S_s S_e S_n S_w real

% Temperatures
syms  T_E T_W T_S T_N T_NW T_NE T_SW T_SE T_P real

% inner Temperatures
syms T_sw T_se T_ne T_nw real

% Define inner Temperatures as interpolation of outer Temperatures
T_sw=(T_SW+T_S+T_P+T_W)/4;
T_se=(T_S+T_SE+T_E+T_P)/4;
T_ne=(T_P+T_E+T_NE+T_N)/4;
T_nw=(T_W+T_P+T_N+T_NW)/4;

% Gradients (Greens theorem)

dTdx_s=  (dy_Sw_Se*T_S + dy_Se_e*T_se + dy_e_w*T_P + dy_w_Sw*T_sw)/S_s;
dTdy_s= -(dx_Sw_Se*T_S + dx_Se_e*T_se + dx_e_w*T_P + dx_w_Sw*T_sw)/S_s;

dTdx_e=  (dy_s_sE*T_se + dy_sE_nE*T_E + dy_nE_n*T_ne + dy_n_s*T_P)/S_e;
dTdy_e= -(dx_s_sE*T_se + dx_sE_nE*T_E + dx_nE_n*T_ne + dx_n_s*T_P)/S_e;

dTdx_n=  (dy_w_e*T_P + dy_e_Ne*T_ne + dy_Ne_Nw*T_N + dy_Nw_w*T_nw)/S_n;
dTdy_n= -(dx_w_e*T_P + dx_e_Ne*T_ne + dx_Ne_Nw*T_N + dx_Nw_w*T_nw)/S_n;

dTdx_w=  (dy_sW_s*T_sw + dy_s_n*T_P + dy_n_nW*T_nw + dy_nW_sW*T_W)/S_w;
dTdy_w= -(dx_sW_s*T_sw + dx_s_n*T_P + dx_n_nW*T_nw + dx_nW_sW*T_W)/S_w;


% Build hole stecil acounting for quadratic lambda like in Helmholtz 

 DDT=(  dy_sw_se*dTdx_s - dx_sw_se*dTdy_s ...
      + dy_se_ne*dTdx_e - dx_se_ne*dTdy_e ...
      + dy_ne_nw*dTdx_n - dx_ne_nw*dTdy_n ...
      + dy_nw_sw*dTdx_w - dx_nw_sw*dTdy_w)/S_P;
 
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
        fprintf(fileID2,'%%    NW(i-1,j-1)   Nw -  N(i-1,j) -  Ne     NE(i-1,j+1)\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%%                 |                 |\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%%       nW - - - - nw ------ n ------ ne - - - nE\n');
        fprintf(fileID2,'%%                 |                 |\n');
        fprintf(fileID2,'%%       |         |        |        |       |\n');
        fprintf(fileID2,'%%                 |                 |\n');
        fprintf(fileID2,'%%   W(i, j-1) - - w - - P (i,j) - - e - -  E (i,j+1)\n');
        fprintf(fileID2,'%%                 |                 |\n');
        fprintf(fileID2,'%%       |         |        |        |       |\n');
        fprintf(fileID2,'%%                 |                 |\n');
        fprintf(fileID2,'%%      sW - - - - sw ------ s ------ se - - - sE\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%%                 |                 |\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%%   SW(i+1,j-1)   Sw  -  S(i+1,j)  - Se      SE(i+1,j+1)\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%% Indexing of stecil: \n\n');
        fprintf(fileID2,'%%    D_4 - D_1 - D2\n');
        fprintf(fileID2,'%%     |     |     | \n\');
        fprintf(fileID2,'%%    D_3 - D_0 - D3\n');
        fprintf(fileID2,'%%     |     |     | \n\');
        fprintf(fileID2,'%%    D_2 -  D1 - D4\n\n');

        
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
