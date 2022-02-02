

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
D3=((dx_e_ne*(dx_P_E/2 + (3*dx_E_nE)/4 + dx_nE_n/4))/S_etae + (dy_e_ne*(dy_P_E/2 + (3*dy_E_nE)/4 + dy_nE_n/4))/S_etae + (dx_e_Ne*dx_ne_nw)/(4*S_n) + (dy_e_Ne*dy_ne_nw)/(4*S_n))/S_eta; 

% West 
D_3=((dx_nw_w*(dx_W_P/2 + (3*dx_nW_W)/4 + dx_n_nW/4))/S_etaw + (dy_nw_w*(dy_W_P/2 + (3*dy_nW_W)/4 + dy_n_nW/4))/S_etaw + (dx_Nw_w*dx_ne_nw)/(4*S_n) + (dy_Nw_w*dy_ne_nw)/(4*S_n))/S_eta; 

% South 
D1=0; 

% North 
D_1=((dx_e_ne*(dx_n_P/4 + dx_nE_n/4))/S_etae + (dx_nw_w*(dx_P_n/4 + dx_n_nW/4))/S_etaw + (dy_e_ne*(dy_n_P/4 + dy_nE_n/4))/S_etae + (dy_nw_w*(dy_P_n/4 + dy_n_nW/4))/S_etaw + (dx_ne_nw*(dx_e_Ne/4 + dx_Nw_w/4 + dx_Ne_Nw))/S_n + (dy_ne_nw*(dy_e_Ne/4 + dy_Nw_w/4 + dy_Ne_Nw))/S_n)/S_eta; 

% NW 
D_4=((dx_nw_w*(dx_nW_W/4 + dx_n_nW/4))/S_etaw + (dy_nw_w*(dy_nW_W/4 + dy_n_nW/4))/S_etaw + (dx_Nw_w*dx_ne_nw)/(4*S_n) + (dy_Nw_w*dy_ne_nw)/(4*S_n))/S_eta; 

% NE 
D2=((dx_e_ne*(dx_E_nE/4 + dx_nE_n/4))/S_etae + (dy_e_ne*(dy_E_nE/4 + dy_nE_n/4))/S_etae + (dx_e_Ne*dx_ne_nw)/(4*S_n) + (dy_e_Ne*dy_ne_nw)/(4*S_n))/S_eta; 

% SW 
D_2=0; 

% SE 
D4=0; 

% P 
D0=((dx_e_ne*(dx_P_E/2 + (3*dx_n_P)/4 + dx_nE_n/4))/S_etae + (dx_nw_w*(dx_W_P/2 + (3*dx_P_n)/4 + dx_n_nW/4))/S_etaw + (dy_e_ne*(dy_P_E/2 + (3*dy_n_P)/4 + dy_nE_n/4))/S_etae + (dy_nw_w*(dy_W_P/2 + (3*dy_P_n)/4 + dy_n_nW/4))/S_etaw + (dx_ne_nw*(dx_w_e + dx_e_Ne/4 + dx_Nw_w/4))/S_n + (dy_ne_nw*(dy_w_e + dy_e_Ne/4 + dy_Nw_w/4))/S_n)/S_eta; 

