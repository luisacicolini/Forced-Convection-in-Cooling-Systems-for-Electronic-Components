
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
target_file = 'build_southEast.m';

fclose(fopen(target_file, 'w'))

%% First part
% Generate stencil with variables names introduced in Camilo F. Silva's
% Course "Numerical Thermo Fluid Dynamics"

% Initialize symbolic variables for distances,
% areas, lambdas, Temperatures

% Around nomega
syms dy_w_P dy_P_N dy_N_Nw dy_Nw_w  real
syms dx_w_P dx_P_N dx_N_Nw dx_Nw_w  real

% Around etaw

syms dy_W_P dy_P_n dy_n_nW dy_nW_W real
syms dx_W_P dx_P_n dx_n_nW dx_nW_W real

% Around P
syms        dy_n_nw dy_nw_w   real
syms        dx_n_nw dx_nw_w   real


% Areas
syms  S_etaomega S_nomega S_etaw   real

% Temperatures
syms  T_E T_W T_S T_N T_NW T_NE T_SW T_SE T_P real

% inner Temperatures
syms T_sw T_se T_ne T_nw real

% Temperatures at the boundaries
syms T_etaW T_etaw T_etaomega T_eta real 
syms T_Nomega T_nomega real

% Define inner Temperatures as interpolation of outer Temperatures
T_sw=(T_SW+T_S+T_P+T_W)/4;
T_se=(T_S+T_SE+T_E+T_P)/4;
T_ne=(T_P+T_E+T_NE+T_N)/4;
T_nw=(T_W+T_P+T_N+T_NW)/4;

T_n     = 0.5*T_P + 0.5*T_N;
T_w     = 0.5*T_P + 0.5*T_W;
T_Nw    = 0.5*T_N + 0.5*T_NW;
T_nW    = 0.5*T_W + 0.5*T_NW;

T_Nomega=   0.75*T_N + 0.25*T_NW;
T_nomega=   0.75*T_n + 0.25*T_nW;
T_omega=    0.75*T_P + 0.25*T_W; 

T_etaW =    0.75*T_W + 0.25*T_NW;
T_etaw =    0.75*T_w + 0.25*T_Nw;
T_etaomega= 0.75*T_omega + 0.25*T_Nomega;
T_eta =     0.75*T_P + 0.25*T_N;


% Gradients (Greens theorem)

dTdx_nomega=  (dy_w_P*T_omega + dy_P_N*T_n + dy_N_Nw*T_Nomega + dy_Nw_w*T_nw)/S_nomega;
dTdy_nomega= -(dx_w_P*T_omega + dx_P_N*T_n + dx_N_Nw*T_Nomega + dx_Nw_w*T_nw)/S_nomega;

dTdx_etaw=  (dy_W_P*T_w + dy_P_n*T_eta + dy_n_nW*T_nw + dy_nW_W*T_etaW)/S_etaw;
dTdy_etaw= -(dx_W_P*T_w + dx_P_n*T_eta + dx_n_nW*T_nw + dx_nW_W*T_etaW)/S_etaw;


% Build hole stecil acounting for quadratic lambda like in Helmholtz 

  DDT =( dy_n_nw*dTdx_nomega - dx_n_nw*dTdy_nomega ...
       + dy_nw_w*dTdx_etaw  - dx_nw_w*dTdy_etaw)/S_etaomega;
 
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
        fprintf(fileID2,'%%    NW=(i-1,j-1)  Nw - N=(i-1,j)- Ne    NE=(i-1,j+1)\n');
        fprintf(fileID2,'%% \n');  
        fprintf(fileID2,'%%                  |                |\n');
        fprintf(fileID2,'%% \n');
        fprintf(fileID2,'%%       nW- - - -  nw----- n ------ne - - - nE\n');
        fprintf(fileID2,'%%                  |                |\n');
        fprintf(fileID2,'%%      etaW        etaw   eta    etae      etaE\n');
        fprintf(fileID2,'%%                  |                |\n');
        fprintf(fileID2,'%%    W=(i,j-1) - - w - - P=(i,j) - -e - - E=(i,j+1)\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%%   SW(i+1,j-1)   Sw  -  S(i+1,j)  - Se      SE(i+1,j+1)\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%% Indexing of stecil: \n\n');
        fprintf(fileID2,'%%    D_4 - D_1 - D2\n');
        fprintf(fileID2,'%%     |     |     | \n\');
        fprintf(fileID2,'%%    D_3 - D_0 - D3\n');

        
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
