

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
D3=((dx_se_e*(dx_E_P/2 + (3*dx_sE_E)/4 + dx_s_sE/4))/S_sigmae + (dy_se_e*(dy_E_P/2 + (3*dy_sE_E)/4 + dy_s_sE/4))/S_sigmae + (dx_s_se*(dx_e_P/4 + dx_Se_e/4))/S_sepsilon + (dy_s_se*(dy_e_P/4 + dy_Se_e/4))/S_sepsilon)/S_sigmaepsilon; 

% West 
D_3=0; 

% South 
D1=((dx_s_se*(dx_P_S/2 + (3*dx_S_Se)/4 + dx_Se_e/4))/S_sepsilon + (dy_s_se*(dy_P_S/2 + (3*dy_S_Se)/4 + dy_Se_e/4))/S_sepsilon + (dx_se_e*(dx_P_s/4 + dx_s_sE/4))/S_sigmae + (dy_se_e*(dy_P_s/4 + dy_s_sE/4))/S_sigmae)/S_sigmaepsilon; 

% North 
D_1=0; 

% NW 
D_4=0; 

% NE 
D2=0; 

% SW 
D_2=0; 

% SE 
D4=((dx_se_e*(dx_sE_E/4 + dx_s_sE/4))/S_sigmae + (dx_s_se*(dx_S_Se/4 + dx_Se_e/4))/S_sepsilon + (dy_se_e*(dy_sE_E/4 + dy_s_sE/4))/S_sigmae + (dy_s_se*(dy_S_Se/4 + dy_Se_e/4))/S_sepsilon)/S_sigmaepsilon; 

% P 
D0=((dx_se_e*(dx_E_P/2 + (3*dx_P_s)/4 + dx_s_sE/4))/S_sigmae + (dx_s_se*(dx_P_S/2 + (3*dx_e_P)/4 + dx_Se_e/4))/S_sepsilon + (dy_se_e*(dy_E_P/2 + (3*dy_P_s)/4 + dy_s_sE/4))/S_sigmae + (dy_s_se*(dy_P_S/2 + (3*dy_e_P)/4 + dy_Se_e/4))/S_sepsilon)/S_sigmaepsilon; 

