

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
D3=0; 

% West 
D_3=((dx_w_sw*(dx_P_W/2 + (3*dx_W_sW)/4 + dx_sW_s/4))/S_sigmaw + (dy_w_sw*(dy_P_W/2 + (3*dy_W_sW)/4 + dy_sW_s/4))/S_sigmaw + (dx_sw_s*(dx_P_w/4 + dx_w_Sw/4))/S_somega + (dy_sw_s*(dy_P_w/4 + dy_w_Sw/4))/S_somega)/S_sigmaomega; 

% South 
D1=((dx_sw_s*(dx_S_P/2 + (3*dx_Sw_S)/4 + dx_w_Sw/4))/S_somega + (dy_sw_s*(dy_S_P/2 + (3*dy_Sw_S)/4 + dy_w_Sw/4))/S_somega + (dx_w_sw*(dx_s_P/4 + dx_sW_s/4))/S_sigmaw + (dy_w_sw*(dy_s_P/4 + dy_sW_s/4))/S_sigmaw)/S_sigmaomega; 

% North 
D_1=0; 

% NW 
D_4=0; 

% NE 
D2=0; 

% SW 
D_2=((dx_sw_s*(dx_Sw_S/4 + dx_w_Sw/4))/S_somega + (dx_w_sw*(dx_W_sW/4 + dx_sW_s/4))/S_sigmaw + (dy_sw_s*(dy_Sw_S/4 + dy_w_Sw/4))/S_somega + (dy_w_sw*(dy_W_sW/4 + dy_sW_s/4))/S_sigmaw)/S_sigmaomega; 

% SE 
D4=0; 

% P 
D0=((dx_sw_s*(dx_S_P/2 + (3*dx_P_w)/4 + dx_w_Sw/4))/S_somega + (dx_w_sw*(dx_P_W/2 + (3*dx_s_P)/4 + dx_sW_s/4))/S_sigmaw + (dy_sw_s*(dy_S_P/2 + (3*dy_P_w)/4 + dy_w_Sw/4))/S_somega + (dy_w_sw*(dy_P_W/2 + (3*dy_s_P)/4 + dy_sW_s/4))/S_sigmaw)/S_sigmaomega; 

