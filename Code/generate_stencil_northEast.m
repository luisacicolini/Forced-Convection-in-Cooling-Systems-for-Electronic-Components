
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
target_file = 'build_northEast.m';

fclose(fopen(target_file, 'w'))

%% First part
% Generate stencil with variables names introduced in Camilo F. Silva's
% Course "Numerical Thermo Fluid Dynamics"

% Initialize symbolic variables for distances,
% areas, lambdas, Temperatures


% Around somega
syms dy_Sw_S dy_S_P dy_P_w dy_w_Sw real
syms dx_Sw_S dx_S_P dx_P_w dx_w_Sw real

% Around sigmaw
syms  dy_sW_s  dy_s_P dy_P_W  dy_W_sW  real
syms  dx_sW_s  dx_s_P dx_P_W  dx_W_sW  real

% Around P

syms dy_sw_s  dy_w_sw     real
syms dx_sw_s  dx_w_sw     real

% Around somega


% Areas
syms  S_sigmaomega S_sigmaw  S_somega   real

% Temperatures
syms  T_E T_W T_S T_N T_NW T_NE T_SW T_SE T_P real

% inner Temperatures
syms T_sw T_se T_ne T_nw real

% Temperatures at the boundaries
syms T_sigmaE T_sigma T_sigmaW real

% Define inner Temperatures as interpolation of outer Temperatures
T_sw=(T_SW+T_S+T_P+T_W)/4;

T_w    = 0.5*T_P + 0.5*T_W;

T_s    = 0.5*T_S + 0.5*T_P;
T_Somega= 0.75*T_S + 0.25*T_SW;
T_omega=  0.75*T_P + 0.25*T_W; 

% T_w    = 0.5*T_P + 0.5*T_W;
T_sigmaW = 0.75*T_W + 0.25*T_SW;
T_sigma = 0.75*T_P + 0.25*T_S;



T_sigmaomega= 0.75*T_omega + 0.25*T_Somega;

% Gradients (Greens theorem)

dTdx_somega=  (dy_Sw_S*T_Somega+dy_S_P*T_s+dy_P_w*T_omega+dy_w_Sw*T_sw)/S_somega;
dTdy_somega= -(dx_Sw_S*T_Somega+dx_S_P*T_s+dx_P_w*T_omega+dx_w_Sw*T_sw)/S_somega;

dTdx_sigmaw=  (dy_sW_s*T_sw+dy_s_P*T_sigma+dy_P_W*T_w+dy_W_sW*T_sigmaW)/S_sigmaw;
dTdy_sigmaw= -(dx_sW_s*T_sw+dx_s_P*T_sigma+dx_P_W*T_w+dx_W_sW*T_sigmaW)/S_sigmaw;


% Build hole stecil acounting for quadratic lambda like in Helmholtz 


  DDT =(  dy_sw_s*dTdx_somega -  dx_sw_s*dTdy_somega ...
       +  dy_w_sw*dTdx_sigmaw -  dx_w_sw*dTdy_sigmaw)/S_sigmaomega;
 
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
