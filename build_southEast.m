

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
D3=0; 

% West 
D_3=((dx_nw_w*(dx_W_P/2 + (3*dx_nW_W)/4 + dx_n_nW/4))/S_etaw + (dy_nw_w*(dy_W_P/2 + (3*dy_nW_W)/4 + dy_n_nW/4))/S_etaw + (dx_n_nw*(dx_w_P/4 + dx_Nw_w/4))/S_nomega + (dy_n_nw*(dy_w_P/4 + dy_Nw_w/4))/S_nomega)/S_etaomega; 

% South 
D1=0; 

% North 
D_1=((dx_n_nw*(dx_P_N/2 + (3*dx_N_Nw)/4 + dx_Nw_w/4))/S_nomega + (dy_n_nw*(dy_P_N/2 + (3*dy_N_Nw)/4 + dy_Nw_w/4))/S_nomega + (dx_nw_w*(dx_P_n/4 + dx_n_nW/4))/S_etaw + (dy_nw_w*(dy_P_n/4 + dy_n_nW/4))/S_etaw)/S_etaomega; 

% NW 
D_4=((dx_nw_w*(dx_nW_W/4 + dx_n_nW/4))/S_etaw + (dx_n_nw*(dx_N_Nw/4 + dx_Nw_w/4))/S_nomega + (dy_nw_w*(dy_nW_W/4 + dy_n_nW/4))/S_etaw + (dy_n_nw*(dy_N_Nw/4 + dy_Nw_w/4))/S_nomega)/S_etaomega; 

% NE 
D2=0; 

% SW 
D_2=0; 

% SE 
D4=0; 

% P 
D0=((dx_nw_w*(dx_W_P/2 + (3*dx_P_n)/4 + dx_n_nW/4))/S_etaw + (dx_n_nw*(dx_P_N/2 + (3*dx_w_P)/4 + dx_Nw_w/4))/S_nomega + (dy_nw_w*(dy_W_P/2 + (3*dy_P_n)/4 + dy_n_nW/4))/S_etaw + (dy_n_nw*(dy_P_N/2 + (3*dy_w_P)/4 + dy_Nw_w/4))/S_nomega)/S_etaomega; 

