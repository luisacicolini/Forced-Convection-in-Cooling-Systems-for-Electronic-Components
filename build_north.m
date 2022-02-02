

% Nomenclature:
%
%    W=(i,j-1) --  w - -P=(i,j) - - e - - E=(i,j+1)
%                  |                |
%   sigmaW    sigmaw    sigma    sigmae    sigmaE
%                  |                |
%      sW - - - -  sw----- s ------se - - - sE
%
%                  |                |
%
%   SW=(i+1,j-1)   Sw  S=(i+1,j) - Se     SE=(i+1,j+1)
%
% Indexing of stecil: 

%    D_3 - D_0 - D3
%     |     |     | 
%    D_2 - D_1 - D4
% Stecil 

% East 
D3=((dx_se_e*(dx_E_P/2 + (3*dx_sE_E)/4 + dx_s_sE/4))/S_sigmae + (dy_se_e*(dy_E_P/2 + (3*dy_sE_E)/4 + dy_s_sE/4))/S_sigmae + (dx_Se_e*dx_sw_se)/(4*S_s) + (dy_Se_e*dy_sw_se)/(4*S_s))/S_sigma; 

% West 
D_3=((dx_w_sw*(dx_P_W/2 + (3*dx_W_sW)/4 + dx_sW_s/4))/S_sigmaw + (dy_w_sw*(dy_P_W/2 + (3*dy_W_sW)/4 + dy_sW_s/4))/S_sigmaw + (dx_w_Sw*dx_sw_se)/(4*S_s) + (dy_w_Sw*dy_sw_se)/(4*S_s))/S_sigma; 

% South 
D1=((dx_se_e*(dx_P_s/4 + dx_s_sE/4))/S_sigmae + (dx_w_sw*(dx_s_P/4 + dx_sW_s/4))/S_sigmaw + (dy_se_e*(dy_P_s/4 + dy_s_sE/4))/S_sigmae + (dy_w_sw*(dy_s_P/4 + dy_sW_s/4))/S_sigmaw + (dx_sw_se*(dx_Se_e/4 + dx_w_Sw/4 + dx_Sw_Se))/S_s + (dy_sw_se*(dy_Se_e/4 + dy_w_Sw/4 + dy_Sw_Se))/S_s)/S_sigma; 

% North 
D_1=0; 

% NW 
D_4=0; 

% NE 
D2=0; 

% SW 
D_2=((dx_w_sw*(dx_W_sW/4 + dx_sW_s/4))/S_sigmaw + (dy_w_sw*(dy_W_sW/4 + dy_sW_s/4))/S_sigmaw + (dx_w_Sw*dx_sw_se)/(4*S_s) + (dy_w_Sw*dy_sw_se)/(4*S_s))/S_sigma; 

% SE 
D4=((dx_se_e*(dx_sE_E/4 + dx_s_sE/4))/S_sigmae + (dy_se_e*(dy_sE_E/4 + dy_s_sE/4))/S_sigmae + (dx_Se_e*dx_sw_se)/(4*S_s) + (dy_Se_e*dy_sw_se)/(4*S_s))/S_sigma; 

% P 
D0=((dx_se_e*(dx_E_P/2 + (3*dx_P_s)/4 + dx_s_sE/4))/S_sigmae + (dx_w_sw*(dx_P_W/2 + (3*dx_s_P)/4 + dx_sW_s/4))/S_sigmaw + (dy_se_e*(dy_E_P/2 + (3*dy_P_s)/4 + dy_s_sE/4))/S_sigmae + (dy_w_sw*(dy_P_W/2 + (3*dy_s_P)/4 + dy_sW_s/4))/S_sigmaw + (dx_sw_se*(dx_e_w + dx_Se_e/4 + dx_w_Sw/4))/S_s + (dy_sw_se*(dy_e_w + dy_Se_e/4 + dy_w_Sw/4))/S_s)/S_sigma; 

