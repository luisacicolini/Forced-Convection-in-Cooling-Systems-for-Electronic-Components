
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
target_file = 'build_southWest.m';

fclose(fopen(target_file, 'w'))

%% First part
% Generate stencil with variables names introduced in Camilo F. Silva's
% Course "Numerical Thermo Fluid Dynamics"

% Initialize symbolic variables for distances,
% areas, lambdas, Temperatures

% Around nepsilon
syms dy_P_e dy_e_Ne dy_Ne_N dy_N_P  real            
syms dx_P_e dx_e_Ne dx_Ne_N dx_N_P  real


% Around etae
syms dy_P_E dy_E_nE dy_nE_n dy_n_P real
syms dx_P_E dx_E_nE dx_nE_n dx_n_P real


% Around P
syms        dy_e_ne dy_ne_n   real
syms        dx_e_ne dx_ne_n   real


% Areas
syms  S_etaepsilon S_nepsilon S_etae   real

% Temperatures
syms  T_E T_W T_S T_N T_NW T_NE T_SW T_SE T_P real

% inner Temperatures
syms T_sw T_se T_ne T_nw real

% Temperatures at the boundaries
syms T_etaE T_etae T_etaepsilon T_eta real 
syms T_Nepsilon T_nepsilon real

% Define inner Temperatures as interpolation of outer Temperatures
T_sw=(T_SW+T_S+T_P+T_W)/4;
T_se=(T_S+T_SE+T_E+T_P)/4;
T_ne=(T_P+T_E+T_NE+T_N)/4;
T_nw=(T_W+T_P+T_N+T_NW)/4;

T_n     = 0.5*T_P + 0.5*T_N;
T_e     = 0.5*T_P + 0.5*T_E;
T_Ne    = 0.5*T_N + 0.5*T_NE;
T_nE    = 0.5*T_E + 0.5*T_NE;

T_Nepsilon=   0.75*T_N + 0.25*T_NE;
T_nepsilon=   0.75*T_n + 0.25*T_nE;
T_epsilon=    0.75*T_P + 0.25*T_E; 

T_etaE =    0.75*T_E + 0.25*T_NE;
T_etae =    0.75*T_e + 0.25*T_Ne;
T_etaepsilon= 0.75*T_epsilon + 0.25*T_Nepsilon;
T_eta =     0.75*T_P + 0.25*T_N;


% Gradients (Greens theorem)

dTdx_nepsilon=  (dy_P_e*T_epsilon + dy_e_Ne*T_ne + dy_Ne_N*T_Nepsilon + dy_N_P*T_n)/S_nepsilon;
dTdy_nepsilon= -(dx_P_e*T_epsilon + dx_e_Ne*T_ne + dx_Ne_N*T_Nepsilon + dx_N_P*T_n)/S_nepsilon;

dTdx_etae=  (dy_P_E*T_e + dy_E_nE*T_etaE + dy_nE_n*T_ne + dy_n_P*T_eta)/S_etae;
dTdy_etae= -(dx_P_E*T_e + dx_E_nE*T_etaE + dx_nE_n*T_ne + dx_n_P*T_eta)/S_etae;


% Build hole stecil acounting for quadratic lambda like in Helmholtz 

  DDT =( dy_e_ne*dTdx_etae - dx_e_ne*dTdy_etae ...
       + dy_ne_n*dTdx_nepsilon  - dx_ne_n*dTdy_nepsilon)/S_etaepsilon;
 
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
