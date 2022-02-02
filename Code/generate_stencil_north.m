
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
target_file = 'build_north.m';

fclose(fopen(target_file, 'w'))

%% First part
% Generate stencil with variables names introduced in Camilo F. Silva's
% Course "Numerical Thermo Fluid Dynamics"

% Initialize symbolic variables for distances,
% areas, lambdas, Temperatures

% Around s
syms dy_Sw_Se dy_Se_e dy_e_w dy_w_Sw real
syms dx_Sw_Se dx_Se_e dx_e_w dx_w_Sw real

% Around sigmae
syms dy_s_sE dy_sE_E dy_E_P dy_P_s  real
syms dx_s_sE dx_sE_E dx_E_P dx_P_s  real

% Around P

syms dy_sw_se  dy_se_e  dy_e_w dy_w_sw     real
syms dx_sw_se  dx_se_e  dx_e_w dx_w_sw     real

% Around sigmaw
syms  dy_sW_s  dy_s_P dy_P_W  dy_W_sW  real
syms  dx_sW_s  dx_s_P dx_P_W  dx_W_sW  real

% Areas
syms S_s S_sigmae S_sigmaw  S_sigma   real

% Temperatures
syms  T_E T_W T_S T_N T_NW T_NE T_SW T_SE T_P real

% inner Temperatures
syms T_sw T_se T_ne T_nw real

% Temperatures at the boundaries
syms T_sigmaE T_sigma T_sigmaW real

% Define inner Temperatures as interpolation of outer Temperatures
T_sw=(T_SW+T_S+T_P+T_W)/4;
T_se=(T_S+T_SE+T_E+T_P)/4;
T_ne=(T_P+T_E+T_NE+T_N)/4;
T_nw=(T_W+T_P+T_N+T_NW)/4;

T_e    = 0.5*T_P + 0.5*T_E;
T_w    = 0.5*T_P + 0.5*T_W;
T_sigmaW = 0.75*T_W + 0.25*T_SW;
T_sigma =  0.75*T_P + 0.25*T_S;
T_sigmaE = 0.75*T_E + 0.25*T_SE;

% Gradients (Greens theorem)

dTdx_s=  (dy_Sw_Se*T_S + dy_Se_e*T_se + dy_e_w*T_P + dy_w_Sw*T_sw)/S_s;
dTdy_s= -(dx_Sw_Se*T_S + dx_Se_e*T_se + dx_e_w*T_P + dx_w_Sw*T_sw)/S_s;

dTdx_sigmae=  (dy_s_sE*T_se + dy_sE_E*T_sigmaE + dy_E_P*T_e + dy_P_s*T_sigma)/S_sigmae;
dTdy_sigmae= -(dx_s_sE*T_se + dx_sE_E*T_sigmaE + dx_E_P*T_e + dx_P_s*T_sigma)/S_sigmae;

dTdx_sigmaw=  (dy_sW_s*T_sw + dy_s_P*T_sigma + dy_P_W*T_w + dy_W_sW*T_sigmaW)/S_sigmaw;
dTdy_sigmaw= -(dx_sW_s*T_sw + dx_s_P*T_sigma + dx_P_W*T_w + dx_W_sW*T_sigmaW)/S_sigmaw;


% Build hole stecil acounting for quadratic lambda like in Helmholtz 


  DDT =( dy_sw_se*dTdx_s     - dx_sw_se*dTdy_s ...
       + dy_se_e*dTdx_sigmae - dx_se_e*dTdy_sigmae ...
       + dy_w_sw*dTdx_sigmaw - dx_w_sw*dTdy_sigmaw)/S_sigma;
 
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
        fprintf(fileID2,'%%    W=(i,j-1) --  w - -P=(i,j) - - e - - E=(i,j+1)\n');
        fprintf(fileID2,'%%                  |                |\n');
        fprintf(fileID2,'%%   sigmaW    sigmaw    sigma    sigmae    sigmaE\n');
        fprintf(fileID2,'%%                  |                |\n');
        fprintf(fileID2,'%%      sW - - - -  sw----- s ------se - - - sE\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%%                  |                |\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%%   SW=(i+1,j-1)   Sw  S=(i+1,j) - Se     SE=(i+1,j+1)\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%% Indexing of stecil: \n\n');
        fprintf(fileID2,'%%    D_3 - D_0 - D3\n');
        fprintf(fileID2,'%%     |     |     | \n\');
        fprintf(fileID2,'%%    D_2 - D_1 - D4\n');

        
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
