

% Nomenclature:
%
%    NW(i-1,j-1)   Nw -  N(i-1,j) -  Ne     NE(i-1,j+1)
%
%                 |                 |
%
%       nW - - - - nw ------ n ------ ne - - - nE
%                 |                 |
%       |         |        |        |       |
%                 |                 |
%   W(i, j-1) - - w - - P (i,j) - - e - -  E (i,j+1)
%                 |                 |
%       |         |        |        |       |
%                 |                 |
%      sW - - - - sw ------ s ------ se - - - sE
%
%                 |                 |
%
%   SW(i+1,j-1)   Sw  -  S(i+1,j)  - Se      SE(i+1,j+1)
%
% Indexing of stecil: 

%    D_4 - D_1 - D2
%     |     |     | 
%    D_3 - D_0 - D3
%     |     |     | 
%    D_2 -  D1 - D4

% Stecil 

% East 
D3=((dx_se_ne*(dx_nE_n/4 + dx_s_sE/4 + dx_sE_nE))/S_e + (dy_se_ne*(dy_nE_n/4 + dy_s_sE/4 + dy_sE_nE))/S_e + (dx_e_Ne*dx_ne_nw)/(4*S_n) + (dx_Se_e*dx_sw_se)/(4*S_s) + (dy_e_Ne*dy_ne_nw)/(4*S_n) + (dy_Se_e*dy_sw_se)/(4*S_s))/S_P; 

% West 
D_3=((dx_nw_sw*(dx_n_nW/4 + dx_sW_s/4 + dx_nW_sW))/S_w + (dy_nw_sw*(dy_n_nW/4 + dy_sW_s/4 + dy_nW_sW))/S_w + (dx_Nw_w*dx_ne_nw)/(4*S_n) + (dx_w_Sw*dx_sw_se)/(4*S_s) + (dy_Nw_w*dy_ne_nw)/(4*S_n) + (dy_w_Sw*dy_sw_se)/(4*S_s))/S_P; 

% South 
D1=((dx_sw_se*(dx_Se_e/4 + dx_w_Sw/4 + dx_Sw_Se))/S_s + (dy_sw_se*(dy_Se_e/4 + dy_w_Sw/4 + dy_Sw_Se))/S_s + (dx_s_sE*dx_se_ne)/(4*S_e) + (dx_sW_s*dx_nw_sw)/(4*S_w) + (dy_s_sE*dy_se_ne)/(4*S_e) + (dy_sW_s*dy_nw_sw)/(4*S_w))/S_P; 

% North 
D_1=((dx_ne_nw*(dx_e_Ne/4 + dx_Nw_w/4 + dx_Ne_Nw))/S_n + (dy_ne_nw*(dy_e_Ne/4 + dy_Nw_w/4 + dy_Ne_Nw))/S_n + (dx_nE_n*dx_se_ne)/(4*S_e) + (dx_n_nW*dx_nw_sw)/(4*S_w) + (dy_nE_n*dy_se_ne)/(4*S_e) + (dy_n_nW*dy_nw_sw)/(4*S_w))/S_P; 

% NW 
D_4=((dx_Nw_w*dx_ne_nw)/(4*S_n) + (dx_n_nW*dx_nw_sw)/(4*S_w) + (dy_Nw_w*dy_ne_nw)/(4*S_n) + (dy_n_nW*dy_nw_sw)/(4*S_w))/S_P; 

% NE 
D2=((dx_nE_n*dx_se_ne)/(4*S_e) + (dx_e_Ne*dx_ne_nw)/(4*S_n) + (dy_nE_n*dy_se_ne)/(4*S_e) + (dy_e_Ne*dy_ne_nw)/(4*S_n))/S_P; 

% SW 
D_2=((dx_w_Sw*dx_sw_se)/(4*S_s) + (dx_sW_s*dx_nw_sw)/(4*S_w) + (dy_w_Sw*dy_sw_se)/(4*S_s) + (dy_sW_s*dy_nw_sw)/(4*S_w))/S_P; 

% SE 
D4=((dx_s_sE*dx_se_ne)/(4*S_e) + (dx_Se_e*dx_sw_se)/(4*S_s) + (dy_s_sE*dy_se_ne)/(4*S_e) + (dy_Se_e*dy_sw_se)/(4*S_s))/S_P; 

% P 
D0=((dx_se_ne*(dx_n_s + dx_nE_n/4 + dx_s_sE/4))/S_e + (dx_ne_nw*(dx_w_e + dx_e_Ne/4 + dx_Nw_w/4))/S_n + (dx_sw_se*(dx_e_w + dx_Se_e/4 + dx_w_Sw/4))/S_s + (dx_nw_sw*(dx_s_n + dx_n_nW/4 + dx_sW_s/4))/S_w + (dy_se_ne*(dy_n_s + dy_nE_n/4 + dy_s_sE/4))/S_e + (dy_ne_nw*(dy_w_e + dy_e_Ne/4 + dy_Nw_w/4))/S_n + (dy_sw_se*(dy_e_w + dy_Se_e/4 + dy_w_Sw/4))/S_s + (dy_nw_sw*(dy_s_n + dy_n_nW/4 + dy_sW_s/4))/S_w)/S_P; 

