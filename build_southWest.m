

% Nomenclature:
%
%    NW=(i-1,j-1)  Nw - N=(i-1,j)- Ne    NE=(i-1,j+1)
% 
%                  |                |
% 
%       nW- - - -  nw----- n ------ne - - - nE
%                  |                |
%      etaW        etaw   eta    etae      etaE
%                  |                |
%    W=(i,j-1) - - w - - P=(i,j) - -e - - E=(i,j+1)
%
%   SW(i+1,j-1)   Sw  -  S(i+1,j)  - Se      SE(i+1,j+1)
%
% Indexing of stecil: 

%    D_4 - D_1 - D2
%     |     |     | 
%    D_3 - D_0 - D3
% Stecil 

% East 
D3=((dx_e_ne*(dx_P_E/2 + (3*dx_E_nE)/4 + dx_nE_n/4))/S_etae + (dy_e_ne*(dy_P_E/2 + (3*dy_E_nE)/4 + dy_nE_n/4))/S_etae + (dx_ne_n*(dx_P_e/4 + dx_e_Ne/4))/S_nepsilon + (dy_ne_n*(dy_P_e/4 + dy_e_Ne/4))/S_nepsilon)/S_etaepsilon; 

% West 
D_3=0; 

% South 
D1=0; 

% North 
D_1=((dx_ne_n*(dx_N_P/2 + (3*dx_Ne_N)/4 + dx_e_Ne/4))/S_nepsilon + (dy_ne_n*(dy_N_P/2 + (3*dy_Ne_N)/4 + dy_e_Ne/4))/S_nepsilon + (dx_e_ne*(dx_n_P/4 + dx_nE_n/4))/S_etae + (dy_e_ne*(dy_n_P/4 + dy_nE_n/4))/S_etae)/S_etaepsilon; 

% NW 
D_4=0; 

% NE 
D2=((dx_e_ne*(dx_E_nE/4 + dx_nE_n/4))/S_etae + (dx_ne_n*(dx_Ne_N/4 + dx_e_Ne/4))/S_nepsilon + (dy_e_ne*(dy_E_nE/4 + dy_nE_n/4))/S_etae + (dy_ne_n*(dy_Ne_N/4 + dy_e_Ne/4))/S_nepsilon)/S_etaepsilon; 

% SW 
D_2=0; 

% SE 
D4=0; 

% P 
D0=((dx_e_ne*(dx_P_E/2 + (3*dx_n_P)/4 + dx_nE_n/4))/S_etae + (dx_ne_n*(dx_N_P/2 + (3*dx_P_e)/4 + dx_e_Ne/4))/S_nepsilon + (dy_e_ne*(dy_P_E/2 + (3*dy_n_P)/4 + dy_nE_n/4))/S_etae + (dy_ne_n*(dy_N_P/2 + (3*dy_P_e)/4 + dy_e_Ne/4))/S_nepsilon)/S_etaepsilon; 

