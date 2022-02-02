

% Nomenclature:
%
%   NW=(i-1,j-1)    Nw -Nomega - N=(i-1,j)
%
%                   |            |
%
%       nW   - - - -nw- nomega - n
%                   |
%       |           |            |
%                   |
%    W=(i,j-1) - -  w - omega   P=(i,j)
%                   |
%       |           |            |
%                   |
%       sW - - -    sw- somega - s
%
%                   |            |
%
%    SW=(i+1,j-1)   Sw -Somega S=(i+1,j)
%
% Indexing of stecil: 

%    D_4 - D_1 
%     |     |     | 
%    D_3 - D_0 
%     |     |     | 
%    D_2 -  D1 

% Stecil 

% East 
D3=0; 

% West 
D_3=((dx_n_nw*(dx_w_P/4 + dx_Nw_w/4))/S_nomega + (dx_sw_s*(dx_P_w/4 + dx_w_Sw/4))/S_somega + (dy_n_nw*(dy_w_P/4 + dy_Nw_w/4))/S_nomega + (dy_sw_s*(dy_P_w/4 + dy_w_Sw/4))/S_somega + (dx_nw_sw*(dx_n_nW/4 + dx_sW_s/4 + dx_nW_sW))/S_w + (dy_nw_sw*(dy_n_nW/4 + dy_sW_s/4 + dy_nW_sW))/S_w)/S_omega; 

% South 
D1=((dx_sw_s*(dx_S_P/2 + (3*dx_Sw_S)/4 + dx_w_Sw/4))/S_somega + (dy_sw_s*(dy_S_P/2 + (3*dy_Sw_S)/4 + dy_w_Sw/4))/S_somega + (dx_sW_s*dx_nw_sw)/(4*S_w) + (dy_sW_s*dy_nw_sw)/(4*S_w))/S_omega; 

% North 
D_1=((dx_n_nw*(dx_P_N/2 + (3*dx_N_Nw)/4 + dx_Nw_w/4))/S_nomega + (dy_n_nw*(dy_P_N/2 + (3*dy_N_Nw)/4 + dy_Nw_w/4))/S_nomega + (dx_n_nW*dx_nw_sw)/(4*S_w) + (dy_n_nW*dy_nw_sw)/(4*S_w))/S_omega; 

% NW 
D_4=((dx_n_nw*(dx_N_Nw/4 + dx_Nw_w/4))/S_nomega + (dy_n_nw*(dy_N_Nw/4 + dy_Nw_w/4))/S_nomega + (dx_n_nW*dx_nw_sw)/(4*S_w) + (dy_n_nW*dy_nw_sw)/(4*S_w))/S_omega; 

% NE 
D2=0; 

% SW 
D_2=((dx_sw_s*(dx_Sw_S/4 + dx_w_Sw/4))/S_somega + (dy_sw_s*(dy_Sw_S/4 + dy_w_Sw/4))/S_somega + (dx_sW_s*dx_nw_sw)/(4*S_w) + (dy_sW_s*dy_nw_sw)/(4*S_w))/S_omega; 

% SE 
D4=0; 

% P 
D0=((dx_n_nw*(dx_P_N/2 + (3*dx_w_P)/4 + dx_Nw_w/4))/S_nomega + (dx_sw_s*(dx_S_P/2 + (3*dx_P_w)/4 + dx_w_Sw/4))/S_somega + (dy_n_nw*(dy_P_N/2 + (3*dy_w_P)/4 + dy_Nw_w/4))/S_nomega + (dy_sw_s*(dy_S_P/2 + (3*dy_P_w)/4 + dy_w_Sw/4))/S_somega + (dx_nw_sw*(dx_s_n + dx_n_nW/4 + dx_sW_s/4))/S_w + (dy_nw_sw*(dy_s_n + dy_n_nW/4 + dy_sW_s/4))/S_w)/S_omega; 

