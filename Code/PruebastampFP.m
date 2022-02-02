
function [stecil,b] = PruebastampFP(i, j, S, Z, TD, alpha, Tinf, boundary, k, beta)
%stecil calculate the linear equation for node (i, j)
%
%  input:
%    i         node number in x direction
%    j         node number in Z direction
%    X         x position of the nodes
%    Z         Z position of the nodes
%    b         right-hand side value for node (i,j)
%    alpha     alpha
%    Tinf      Tinf for Robin BC
%    boundary  defines the boundary conditions
%    verbose   verbositiy level
%
%  output:
%    stecil     linear equation for node (i,j)
%    b         new right-hand side value for node (i,j)


% Init

b=0;
n = size(S, 1);
m = size(S, 2);
stecil = zeros(1, n*m);
index=@(ii, jj) ii + (jj-1)*n;

% Determine the node positon
nodePosition = 'inner Node';


if i==n
    nodePosition = 'South';
elseif i==1
    nodePosition = 'North';
end


if j == 1
    nodePosition = 'West';
    
    if i == 1
        nodePosition = 'WestNorth';
    elseif i == n
        nodePosition = 'WestSouth';
    end
    
elseif j == m
    nodePosition = 'East';
    
    if i == 1
        nodePosition = 'EastNorth';
    elseif i == n
        nodePosition = 'EastSouth';
    end
end




% Calculate the equation for the correct node position
switch nodePosition
    
    case 'inner Node'
        
        
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
        
        
        
        y_NW = Z(i-1,j-1);   x_NW = S(i-1,j-1);
        y_N  = Z(i-1,j  );   x_N  = S(i-1,j  );
        y_NE = Z(i-1,j+1);   x_NE = S(i-1,j+1);
        y_W  = Z(i  ,j-1);   x_W  = S(i  ,j-1);
        y_P  = Z(i  ,j  );   x_P  = S(i  ,j  );
        y_E  = Z(i  ,j+1);   x_E  = S(i  ,j+1);
        y_SW = Z(i+1,j-1);   x_SW = S(i+1,j-1);
        y_S  = Z(i+1,j  );   x_S  = S(i+1,j  );
        y_SE = Z(i+1,j+1);   x_SE = S(i+1,j+1);
        
        y_Nw = (y_N + y_NW)/2;  x_Nw = (x_N + x_NW)/2;
        y_Ne = (y_NE + y_N)/2;  x_Ne = (x_NE + x_N)/2;
        y_w  = (y_P  + y_W)/2;  x_w  = (x_P  + x_W)/2;
        y_e  = (y_E  + y_P)/2;  x_e  = (x_E  + x_P)/2;
        y_Sw = (y_S  + y_SW)/2; x_Sw = (x_S  + x_SW)/2;
        y_Se = (y_SE + y_S)/2;  x_Se = (x_SE + x_S)/2;
        
        y_nW = (y_NW + y_W)/2;   x_nW = (x_NW + x_W)/2;
        y_sW = (y_W  + y_SW)/2;  x_sW = (x_W  + x_SW)/2;
        y_n  =  (y_N + y_P)/2;   x_n  =  (x_N + x_P)/2;
        y_s  =  (y_P +  y_S)/2;  x_s  =  (x_P +  x_S)/2;
        y_nE =  (y_NE + y_E)/2;  x_nE =  (x_NE + x_E)/2;
        y_sE =  (y_E +  y_SE)/2; x_sE =  (x_E +  x_SE)/2;
        
        y_se = (y_Se + y_e)/2;  x_se = (x_Se + x_e)/2;
        y_sw = (y_w + y_Sw)/2;  x_sw = (x_w + x_Sw)/2;
        y_nw = (y_w + y_Nw)/2;  x_nw = (x_w + x_Nw)/2;
        y_ne = (y_Ne + y_e)/2;  x_ne = (x_Ne + x_e)/2;
        
        
        
        % Around s
        
        dy_Sw_Se = y_Se - y_Sw;  dx_Sw_Se = x_Se - x_Sw;
        dy_Se_e = y_e - y_Se;    dx_Se_e = x_e - x_Se;
        dy_e_w = y_w - y_e;      dx_e_w = x_w - x_e;
        dy_w_Sw = y_Sw - y_w;    dx_w_Sw = x_Sw - x_w;
        
        % Around e
        
        dy_s_sE = y_sE - y_s;    dx_s_sE = x_sE - x_s;
        dy_sE_nE = y_nE - y_sE;  dx_sE_nE = x_nE - x_sE;
        dy_nE_n = y_n - y_nE;    dx_nE_n = x_n - x_nE;
        dy_n_s = y_s - y_n;      dx_n_s = x_s - x_n;
        
        % Around n
        
        dy_w_e   = y_e - y_w;    dx_w_e   = x_e - x_w;
        dy_e_Ne  = y_Ne - y_e;   dx_e_Ne  = x_Ne - x_e;
        dy_Ne_Nw = y_Nw - y_Ne;  dx_Ne_Nw = x_Nw - x_Ne;
        dy_Nw_w  = y_w - y_Nw;   dx_Nw_w  = x_w - x_Nw;
        
        
        % Around w
        dy_sW_s = y_s - y_sW;    dx_sW_s = x_s - x_sW;
        dy_s_n  = y_n - y_s;     dx_s_n  = x_n - x_s;
        dy_n_nW = y_nW - y_n;    dx_n_nW = x_nW - x_n;
        dy_nW_sW = y_sW - y_nW;  dx_nW_sW = x_sW - x_nW;
        
        
        % Around P
        
        dy_sw_se = y_se - y_sw;   dx_sw_se = x_se - x_sw;
        dy_se_ne = y_ne - y_se;   dx_se_ne = x_ne - x_se;
        dy_ne_nw = y_nw - y_ne;   dx_ne_nw = x_nw - x_ne;
        dy_nw_sw = y_sw - y_nw;   dx_nw_sw = x_sw - x_nw;
        
        % Areas
        
        
        
        S_P                 = abs( x_ne*y_se - y_ne*x_se...
            + x_se*y_sw - y_se*x_sw...
            + x_sw*y_nw - y_sw*x_nw...
            + x_nw*y_ne - y_nw*x_ne) / 2;
        S_s                 = abs( x_e*y_Se - y_e*x_Se...
            + x_Se*y_Sw - y_Se*x_Sw...
            + x_Sw*y_w - y_Sw*x_w...
            + x_w*y_e - y_w*x_e) / 2;
        S_e                 = abs( x_nE*y_sE - y_nE*x_sE...
            + x_sE*y_s - y_sE*x_s...
            + x_s*y_n - y_s*x_n...
            + x_n*y_nE - y_n*x_nE) / 2;
        S_n                 = abs( x_Ne*y_e - y_Ne*x_e...
            + x_e*y_w - y_e*x_w...
            + x_w*y_Nw - y_w*x_Nw...
            + x_Nw*y_Ne - y_Nw*x_Ne ) / 2;
        S_w                 = abs( x_n*y_s - y_n*x_s...
            + x_s*y_sW - y_s*x_sW...
            + x_sW*y_nW - y_sW*x_nW...
            + x_nW*y_n - y_nW*x_n) / 2;
        
        
        
        %$$$$$$$$$$$$$$$$$$$$$$ Stecil $$$$$$$$$$$$$$$$$$$
        
        build_inner
        
        %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        
        % East
        stecil(index(i, j+1)) = D3;
        % West
        stecil(index(i, j-1)) = D_3;
        % South
        stecil(index(i+1, j)) = D1;
        % North
        stecil(index(i-1, j)) = D_1;
        % NW
        stecil(index(i-1, j-1)) = D_4;
        % NE
        stecil(index(i-1, j+1)) = D2;
        % SW
        stecil(index(i+1, j-1)) = D_2;
        % SE
        stecil(index(i+1, j+1)) = D4;
        % P
        stecil(index(i, j)) = D0;
        
    case 'South'
        
        if strcmp(boundary.south, 'Neumann')
            
            % Nomencature:
            %    NW=(i-1,j-1)  Nw - N=(i-1,j)- Ne    NE=(i-1,j+1)
            %
            %                  |                |
            %
            %       nW- - - -  nw----- n ------ne - - - nE
            %                  |                |
            %      etaW        etaw   eta    etae      etaE
            %                  |                |
            %    W=(i,j-1) - - w - - P=(i,j) - -e - - E=(i,j+1)
            
            
            
            
            y_NW = Z(i-1,j-1);   x_NW = S(i-1,j-1);
            y_N  = Z(i-1,j  );   x_N  = S(i-1,j  );
            y_NE = Z(i-1,j+1);   x_NE = S(i-1,j+1);
            y_W  = Z(i  ,j-1);   x_W  = S(i  ,j-1);
            y_P  = Z(i  ,j  );   x_P  = S(i  ,j  );
            y_E  = Z(i  ,j+1);   x_E  = S(i  ,j+1);
            
            
            y_Nw = (y_N + y_NW)/2;  x_Nw = (x_N + x_NW)/2;
            y_Ne = (y_NE + y_N)/2;  x_Ne = (x_NE + x_N)/2;
            y_w  = (y_P  + y_W)/2;  x_w  = (x_P  + x_W)/2;
            y_e  = (y_E  + y_P)/2;  x_e  = (x_E  + x_P)/2;
            
            
            y_nW = (y_NW + y_W)/2;   x_nW = (x_NW + x_W)/2;
            y_n  =  (y_N + y_P)/2;   x_n  =  (x_N + x_P)/2;
            y_nE =  (y_NE + y_E)/2;  x_nE =  (x_NE + x_E)/2;
            
            
            y_nw = (y_w + y_Nw)/2;  x_nw = (x_w + x_Nw)/2;
            y_ne = (y_Ne + y_e)/2;  x_ne = (x_Ne + x_e)/2;
            
            y_etaW = (y_nW + y_W)/2; x_etaW = (x_nW + x_W)/2;
            y_etaw = (y_nw + y_w)/2; x_etaw = (x_nw + x_w)/2;
            y_eta  = (y_n  + y_P)/2; x_eta  = (x_n  + x_P)/2;
            y_etae = (y_ne + y_e)/2; x_etae = (x_ne + x_e)/2;
            y_etaE = (y_nE + y_E)/2; x_etaE = (x_nE + x_E)/2;
            
            
            % Around etaw
            
            dy_W_P  = y_P - y_W;    dx_W_P  = x_P - x_W;
            dy_P_n  = y_n - y_P;    dx_P_n  = x_n - x_P;
            dy_n_nW = y_nW - y_n;   dx_n_nW = x_nW - x_n;
            dy_nW_W = y_W - y_nW;   dx_nW_W = x_W - x_nW;
            
            
            % Around etae
            
            dy_P_E  = y_E - y_P;    dx_P_E  = x_E - x_P;
            dy_E_nE = y_nE - y_E;   dx_E_nE = x_nE - x_E;
            dy_nE_n = y_n - y_nE;   dx_nE_n = x_n - x_nE;
            dy_n_P  = y_P - y_n;    dx_n_P  = x_P - x_n;
            
            % Around n
            
            dy_w_e   = y_e - y_w;    dx_w_e   = x_e - x_w;
            dy_e_Ne  = y_Ne - y_e;   dx_e_Ne  = x_Ne - x_e;
            dy_Ne_Nw = y_Nw - y_Ne;  dx_Ne_Nw = x_Nw - x_Ne;
            dy_Nw_w  = y_w - y_Nw;   dx_Nw_w  = x_w - x_Nw;
            
            
            % Around P
            
            dy_e_ne = y_ne - y_e;     dx_e_ne = x_ne - x_e;
            dy_ne_nw = y_nw - y_ne;   dx_ne_nw = x_nw - x_ne;
            dy_nw_w =  y_w - y_nw;    dx_nw_w =  x_w - x_nw;
            
            % At boundary
            dl_w_e = norm([dx_w_e; dy_w_e]);
            
            % Areas
            
            
            S_eta                 = abs( x_ne*y_e - y_ne*x_e...
                + x_e*y_w - y_e*x_w...
                + x_w*y_nw - y_w*x_nw...
                + x_nw*y_ne - y_nw*x_ne) / 2;
            S_etae                 = abs( x_nE*y_E - y_nE*x_E...
                + x_E*y_P - y_E*x_P...
                + x_P*y_n - y_P*x_n...
                + x_n*y_nE - y_n*x_nE) / 2;
            S_n                 = abs( x_Ne*y_e - y_Ne*x_e...
                + x_e*y_w - y_e*x_w...
                + x_w*y_Nw - y_w*x_Nw...
                + x_Nw*y_Ne - y_Nw*x_Ne ) / 2;
            S_etaw                 = abs( x_n*y_P - y_n*x_P...
                + x_P*y_W - y_P*x_W...
                + x_W*y_nW - y_W*x_nW...
                + x_nW*y_n - y_nW*x_n) / 2;
            
            
            %$$$$$$$$$$$$$$$$$$$$$$ Stecil $$$$$$$$$$$$$$$$$$$
            
            build_south
            
            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            
            % East
            stecil(index(i, j+1)) = D3;
            % West
            stecil(index(i, j-1)) = D_3;
            % North
            stecil(index(i-1, j)) = D_1;
            % NW
            stecil(index(i-1, j-1)) = D_4;
            % NE
            stecil(index(i-1, j+1)) = D2;
            % P
            stecil(index(i, j)) = D0;
            
            b = -dl_w_e/S_eta*(-beta);
            
        elseif strcmp(boundary.south, 'Dirichlet')
            stecil(index(i, j)) = 1;
            b=TD.south;
        else
            warning('boundary condition not implemented')
        end
        
    case 'North'
        
        if strcmp(boundary.north, 'Robin') || strcmp(boundary.north, 'Neumann')
            

            
            % Nomencature:
            %    W=(i,j-1) --  w - -P=(i,j) - - e - - E=(i,j+1)
            %                  |                |
            %   sigmaW    sigmaw    sigma    sigmae    sigmaE
            %                  |                |
            %      sW - - - -  sw----- s ------se - - - sE
            %
            %                  |                |
            %
            %   SW=(i+1,j-1)   Sw  S=(i+1,j) - Se     SE=(i+1,j+1)
            
            
            % Calculate some help values
            
            
            y_W  = Z(i  ,j-1);   x_W  = S(i  ,j-1);
            y_P  = Z(i  ,j  );   x_P  = S(i  ,j  );
            y_E  = Z(i  ,j+1);   x_E  = S(i  ,j+1);
            y_SW = Z(i+1,j-1);   x_SW = S(i+1,j-1);
            y_S  = Z(i+1,j  );   x_S  = S(i+1,j  );
            y_SE = Z(i+1,j+1);   x_SE = S(i+1,j+1);
            
            y_w  = (y_P  + y_W)/2;  x_w  = (x_P  + x_W)/2;
            y_e  = (y_E  + y_P)/2;  x_e  = (x_E  + x_P)/2;
            y_Sw = (y_S  + y_SW)/2; x_Sw = (x_S  + x_SW)/2;
            y_Se = (y_SE + y_S)/2;  x_Se = (x_SE + x_S)/2;
            
            y_sW = (y_W  + y_SW)/2;  x_sW = (x_W  + x_SW)/2;
            y_s  =  (y_P +  y_S)/2;  x_s  =  (x_P +  x_S)/2;
            y_sE =  (y_E +  y_SE)/2; x_sE =  (x_E +  x_SE)/2;
            
            y_se = (y_Se + y_e)/2;  x_se = (x_Se + x_e)/2;
            y_sw = (y_w + y_Sw)/2;  x_sw = (x_w + x_Sw)/2;
            
            
            y_sigmaW = (y_W + y_sW)/2; x_sigmaW = (y_W + y_sW)/2;
            y_sigmaw = (y_w + y_sw)/2; x_sigmaw = (y_w + y_sw)/2;
            y_sigma  = (y_P + y_s)/2;  x_sigma  = (y_P + y_s)/2;
            y_sigmae = (y_e + y_se)/2; x_sigmae = (y_e + y_se)/2;
            y_sigmaE = (y_E + y_sE)/2; x_sigmaE = (y_E + y_sE)/2;
            
            
            
            % Around s
            
            dy_Sw_Se = y_Se - y_Sw; dx_Sw_Se = x_Se - x_Sw;
            dy_Se_e  = y_e - y_Se ; dx_Se_e  = x_e - x_Se ;
            dy_e_w   = y_w - y_e  ; dx_e_w   = x_w - x_e  ;
            dy_w_Sw  = y_Sw - y_w ; dx_w_Sw  = x_Sw - x_w ;
            
            % Around sigmae
            
            dy_s_sE  = y_sE - y_s ; dx_s_sE  = x_sE - x_s ;
            dy_sE_E  = y_E - y_sE ; dx_sE_E  = x_E - x_sE ;
            dy_E_P   = y_P - y_E ;  dx_E_P   = x_P - x_E ;
            dy_P_s   = y_s - y_P ;  dx_P_s   = x_s - x_P ;
            
            % Around sigmaw
            
            dy_sW_s  = y_s - y_sW;  dx_sW_s  = x_s - x_sW;
            dy_s_P   = y_P - y_s ;  dx_s_P   = x_P - x_s ;
            dy_P_W   = y_W - y_P ;  dx_P_W   = x_W - x_P ;
            dy_W_sW  = y_sW - y_W;  dx_W_sW  = x_sW - x_W;
            
            % Around P
            
            dy_sw_se  = y_se - y_sw;  dx_sw_se  = x_se - x_sw;
            dy_se_e   = y_e - y_se ;  dx_se_e   = x_e - x_se ;
            dy_e_w    = y_w - y_e;    dx_e_w    = x_w - x_e;
            dy_w_sw   = y_sw - y_w;   dx_w_sw   = x_sw - x_w;
            
            
            % At boundary
            dl_e_w = norm([dx_e_w; dy_e_w]);
            
            
            
            % Areas
            
            S_sigma                 = abs( x_e*y_se - y_e*x_se...
                + x_se*y_sw - y_se*x_sw...
                + x_sw*y_w - y_sw*x_w...
                + x_w*y_e - y_w*x_e) / 2;
            S_s                 = abs( x_e*y_Se - y_e*x_Se...
                + x_Se*y_Sw - y_Se*x_Sw...
                + x_Sw*y_w - y_Sw*x_w...
                + x_w*y_e - y_w*x_e) / 2;
            S_sigmae                 = abs( x_E*y_sE - y_E*x_sE...
                + x_sE*y_s - y_sE*x_s...
                + x_s*y_P - y_s*x_P...
                + x_P*y_E - y_P*x_E) / 2;
            S_sigmaw                 = abs( x_P*y_s - y_P*x_s...
                + x_s*y_sW - y_s*x_sW...
                + x_sW*y_W - y_sW*x_W...
                + x_W*y_P - y_W*x_P) / 2;
            
            %$$$$$$$$$$$$$$$$$$$$$$ Stecil $$$$$$$$$$$$$$$$$$$
            
            build_north
            
            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            
            % East
            stecil(index(i, j+1)) = D3;
            % West
            stecil(index(i, j-1)) = D_3;
            % South
            stecil(index(i+1, j)) = D1;
            % SW
            stecil(index(i+1, j-1)) = D_2;
            % SE
            stecil(index(i+1, j+1)) = D4;
            % P
            stecil(index(i, j)) = D0;
            
            if strcmp(boundary.north, 'Neumann')
                b = -dl_e_w/S_sigma*(-beta);
            elseif strcmp(boundary.north, 'Robin')                
                stecil(index(i,j)) = stecil(index(i, j)) + ...
                    dl_e_w/S_sigma*(-1/k);
                b = -dl_e_w/S_sigma*(-1/k*(-Tinf));
            end
            
        elseif strcmp(boundary.north, 'Dirichlet')
            stecil(index(i,j)) = 1;
            b=TD.north;
        else
            warning('boundary condition not implemented')
        end
        
        
    case 'East'
        
        if strcmp(boundary.east, 'Robin') || strcmp(boundary.east, 'Neumann')
            

            
            % Nomencature:
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
            
            % Calculate some help values
            
            
            y_NW = Z(i-1,j-1);   x_NW = S(i-1,j-1);
            y_N  = Z(i-1,j  );   x_N  = S(i-1,j  );
            y_W  = Z(i  ,j-1);   x_W  = S(i  ,j-1);
            y_P  = Z(i  ,j  );   x_P  = S(i  ,j  );
            y_SW = Z(i+1,j-1);   x_SW = S(i+1,j-1);
            y_S  = Z(i+1,j  );   x_S  = S(i+1,j  );
            
            y_Nw = (y_N + y_NW)/2;  x_Nw = (x_N + x_NW)/2;
            y_w  = (y_P  + y_W)/2;  x_w  = (x_P  + x_W)/2;
            y_Sw = (y_S  + y_SW)/2; x_Sw = (x_S  + x_SW)/2;
            
            y_nW = (y_NW + y_W)/2;   x_nW = (x_NW + x_W)/2;
            y_sW = (y_W  + y_SW)/2;  x_sW = (x_W  + x_SW)/2;
            y_n  =  (y_N + y_P)/2;   x_n  =  (x_N + x_P)/2;
            y_s  =  (y_P +  y_S)/2;  x_s  =  (x_P +  x_S)/2;
            
            y_sw = (y_w + y_Sw)/2;  x_sw = (x_w + x_Sw)/2;
            y_nw = (y_w + y_Nw)/2;  x_nw = (x_w + x_Nw)/2;
            
            
            y_Nomega  = (y_N + y_Nw)/2;  x_Nomega  = (x_N + x_Nw)/2;
            y_nomega  = (y_n  + y_nw)/2; x_nomega  = (x_n  + x_nw)/2;
            y_omega   = (y_P + y_w)/2;   x_omega   = (x_P + x_w)/2;
            y_somega  = (y_s + y_sw)/2;  x_somega  = (x_s + x_sw)/2;
            y_Somega  = (y_S + y_Sw)/2;  x_Somega  = (x_S + x_Sw)/2;
            
            
            % Around somega
            
            dy_Sw_S = y_S - y_Sw;  dx_Sw_S = x_S - x_Sw;
            dy_S_P  = y_P - y_S;   dx_S_P  = x_P - x_S;
            dy_P_w  = y_w - y_P;   dx_P_w  = x_w - x_P;
            dy_w_Sw = y_Sw - y_w;  dx_w_Sw = x_Sw - x_w;
            
            % Around nomega
            
            dy_w_P  = y_P - y_w;   dx_w_P  = x_P - x_w;
            dy_P_N  = y_N - y_P;   dx_P_N  = x_N - x_P;
            dy_N_Nw = y_Nw - y_N;  dx_N_Nw = x_Nw - x_N;
            dy_Nw_w = y_w - y_Nw;  dx_Nw_w = x_w - x_Nw;
            
            % Around w
            
            dy_sW_s  = y_s - y_sW;   dx_sW_s  = x_s - x_sW;
            dy_s_n   = y_n - y_s;    dx_s_n   = x_n - x_s;
            dy_n_nW  = y_nW - y_n;   dx_n_nW  = x_nW - x_n;
            dy_nW_sW = y_sW - y_nW;  dx_nW_sW = x_sW - x_nW;
            
            % Around P
            
            dy_sw_s  = y_s - y_sw;   dx_sw_s  = x_s - x_sw;
            dy_s_n   = y_n - y_s;    dx_s_n   = x_n - x_s;
            dy_n_nw  = y_nw - y_n;   dx_n_nw  = x_nw - x_n;
            dy_nw_sw = y_sw - y_nw;  dx_nw_sw = x_sw - x_nw;
            
            % At boundary
            
            dl_s_n = norm([dx_s_n; dy_s_n]);
            
            % Areas
            
            
            S_omega                 = abs( x_n*y_s - y_n*x_s...
                + x_s*y_sw - y_s*x_sw...
                + x_sw*y_nw - y_sw*x_nw...
                + x_nw*y_n - y_nw*x_n) / 2;
            S_somega                 = abs( x_P*y_S - y_P*x_S...
                + x_S*y_Sw - y_S*x_Sw...
                + x_Sw*y_w - y_Sw*x_w...
                + x_w*y_P - y_w*x_P) / 2;
            
            S_nomega                 = abs( x_N*y_P - y_N*x_P...
                + x_P*y_w - y_P*x_w...
                + x_w*y_Nw - y_w*x_Nw...
                + x_Nw*y_N - y_Nw*x_N ) / 2;
            S_w                 = abs( x_n*y_s - y_n*x_s...
                + x_s*y_sW - y_s*x_sW...
                + x_sW*y_nW - y_sW*x_nW...
                + x_nW*y_n - y_nW*x_n) / 2;
            
            
            %$$$$$$$$$$$$$$$$$$$$$$ Stecil $$$$$$$$$$$$$$$$$$$
            
            build_east
            
            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            
            
            % West
            stecil(index(i, j-1)) = D_3;
            % South
            stecil(index(i+1, j)) = D1;
            % North
            stecil(index(i-1, j)) = D_1;
            % NW
            stecil(index(i-1, j-1)) = D_4;
            % SW
            stecil(index(i+1, j-1)) = D_2;
            % P
            stecil(index(i, j)) = D0 ;
            
            if strcmp(boundary.east, 'Neumann')
                b = -dl_s_n/S_omega*(-beta);
            elseif strcmp(boundary.east, 'Robin')
                stecil(index(i,j)) = D0 + ...
                    dl_s_n/S_omega*(-alpha/k);
                b = -dl_s_n/S_omega*(-alpha/k*(-Tinf));
            end           
            
        elseif strcmp(boundary.east, 'Dirichlet')
            stecil(index(i,j)) = 1;
            b=TD.east;
        else
            warning('boundary condition not implemented')
        end
        
        
    case 'West'
        
        if strcmp(boundary.west, 'Robin') || strcmp(boundary.west, 'Neumann')
            

            
            % Calculate some help values
            
            
            y_NE = Z(i-1,j+1);   x_NE = S(i-1,j+1);
            y_N  = Z(i-1,j  );   x_N  = S(i-1,j  );
            y_E  = Z(i  ,j+1);   x_E  = S(i  ,j+1);
            y_P  = Z(i  ,j  );   x_P  = S(i  ,j  );
            y_SE = Z(i+1,j+1);   x_SE = S(i+1,j+1);
            y_S  = Z(i+1,j  );   x_S  = S(i+1,j  );
            
            y_Ne = (y_N + y_NE)/2;  x_Ne = (x_N + x_NE)/2;
            y_e  = (y_P  + y_E)/2;  x_e  = (x_P  + x_E)/2;
            y_Se = (y_S  + y_SE)/2; x_Se = (x_S  + x_SE)/2;
            
            y_nE = (y_NE + y_E)/2;   x_nE = (x_NE + x_E)/2;
            y_sE = (y_E  + y_SE)/2;  x_sE = (x_E  + x_SE)/2;
            y_n  =  (y_N + y_P)/2;   x_n  =  (x_N + x_P)/2;
            y_s  =  (y_P +  y_S)/2;  x_s  =  (x_P +  x_S)/2;
            
            y_se = (y_e + y_Se)/2;  x_se = (x_e + x_Se)/2;
            y_ne = (y_e + y_Ne)/2;  x_ne = (x_e + x_Ne)/2;
            
            
            y_Nepsilon  = (y_N + y_Ne)/2;  x_Nepsilon  = (x_N + x_Ne)/2;
            y_nepsilon  = (y_n  + y_ne)/2; x_nepsilon  = (x_n  + x_ne)/2;
            y_epsilon   = (y_P + y_e)/2;   x_epsilon   = (x_P + x_e)/2;
            y_sepsilon  = (y_s + y_se)/2;  x_sepsilon  = (x_s + x_se)/2;
            y_Sepsilon  = (y_S + y_Se)/2;  x_Sepsilon  = (x_S + x_Se)/2;
            
            
            % Around sepsilon
            
            dy_S_Se  = y_Se - y_S;   dx_S_Se = x_Se - x_S;
            dy_Se_e  = y_e - y_Se;   dx_Se_e = x_e - x_Se;
            dy_e_P   = y_P - y_e;    dx_e_P  = x_P - x_e;
            dy_P_S   = y_S - y_P;    dx_P_S  = x_S - x_P;
            
            % Around nepsilon
            
            dy_P_e   = y_e - y_P;    dx_P_e   = x_e - x_P;
            dy_e_Ne  = y_Ne - y_e;   dx_e_Ne  = x_Ne - x_e;
            dy_Ne_N  = y_N - y_Ne;   dx_Ne_N  = x_N - x_Ne;
            dy_N_P   = y_P - y_N;    dx_N_P   = x_P - x_N;
            
            % Around e
            
            dy_s_sE  = y_sE - y_s;   dx_s_sE  = x_sE - x_s;
            dy_sE_nE = y_nE - y_sE;  dx_sE_nE = x_nE - x_sE;
            dy_nE_n  = y_n - y_nE;   dx_nE_n  = x_n - x_nE;
            dy_n_s   = y_s - y_n;    dx_n_s   = x_s - x_n;
            
            % Around P
            
            dy_s_se  = y_se - y_s;   dx_s_se  = x_se - x_s;
            dy_se_ne = y_ne - y_se;  dx_se_ne = x_ne - x_se;
            dy_ne_n  = y_n - y_ne;   dx_ne_n  = x_n - x_ne;
            dy_n_s   = y_s - y_n;    dx_n_s   = x_s - x_n;
            
            % At boundary
            
            dl_n_s = norm([dx_n_s; dy_n_s]);
            
            % Areas
            
            
            S_epsilon                 = abs( x_s*y_n - y_s*x_n...
                + x_n*y_ne - y_n*x_ne...
                + x_ne*y_se - y_ne*x_se...
                + x_se*y_s - y_se*x_s) / 2;
            
            S_sepsilon                 = abs( x_S*y_P - y_S*x_P...
                + x_P*y_e - y_P*x_e...
                + x_e*y_Se - y_e*x_Se...
                + x_Se*y_S - y_Se*x_S) / 2;
            
            S_nepsilon                 = abs( x_P*y_N - y_P*x_N...
                + x_N*y_Ne - y_N*x_Ne...
                + x_Ne*y_e - y_Ne*x_e...
                + x_e*y_P - y_e*x_P ) / 2;
            
            S_e                 = abs( x_s*y_n - y_s*x_n...
                + x_n*y_nE - y_n*x_nE...
                + x_nE*y_sE - y_nE*x_sE...
                + x_sE*y_s - y_sE*x_s) / 2;
            
            
            %$$$$$$$$$$$$$$$$$$$$$$ Stecil $$$$$$$$$$$$$$$$$$$
            
            build_west
            
            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            
            % East
            stecil(index(i, j+1)) = D3;
            % South
            stecil(index(i+1, j)) = D1;
            % North
            stecil(index(i-1, j)) = D_1;
            % NE
            stecil(index(i-1, j+1)) = D2;
            % SE
            stecil(index(i+1, j+1)) = D4;
            % P
            stecil(index(i, j)) = D0 ;
            
            
            if strcmp(boundary.west, 'Neumann')
                b = -dl_n_s/S_epsilon*(-beta);
            elseif strcmp(boundary.west, 'Robin')
                stecil(index(i,j)) = D0 + ...
                    dl_n_s/S_epsilon*(-alpha/k);
                b = -dl_n_s/S_epsilon*(-alpha/k*(-Tinf));
            end
            
            
        elseif strcmp(boundary.west, 'Dirichlet')
            stecil(index(i,j)) = 1;
            b=TD.west;
        else
            warning('boundary condition not implemented')
        end
        
        
    case 'EastSouth'
        
        if strcmp(boundary.south, 'Neumann')
            
            if strcmp(boundary.east, 'Neumann') || strcmp(boundary.east, 'Robin')
                
                % Nomencature:
                %    NW=(i-1,j-1)  Nw - Nomega  - N=(i-1,j)
                %                  |       |        |
                %                  |       |        |
                %                  |       |        |
                %       nW- - - -  nw----nomega  ----- n
                %                  |       |        |
                %      etaW      etaw   etaomega   eta
                %                  |       |        |
                %    W=(i,j-1) - - w - - omega -- P=(i,j)
                
                
                % Calculate some help values
                
                y_NW = Z(i-1,j-1);   x_NW = S(i-1,j-1);
                y_N  = Z(i-1,j  );   x_N  = S(i-1,j  );
                y_W  = Z(i  ,j-1);   x_W  = S(i  ,j-1);
                y_P  = Z(i  ,j  );   x_P  = S(i  ,j  );
                
                y_Nw = (y_N + y_NW)/2;  x_Nw = (x_N + x_NW)/2;
                y_w  = (y_P  + y_W)/2;  x_w  = (x_P  + x_W)/2;
                
                y_nW = (y_NW + y_W)/2;   x_nW = (x_NW + x_W)/2;
                y_n  =  (y_N + y_P)/2;   x_n  =  (x_N + x_P)/2;
                
                y_nw = (y_w + y_Nw)/2;  x_nw = (x_w + x_Nw)/2;
                
                % Around nomega
                
                dy_w_P  = y_P - y_w;   dx_w_P  = x_P - x_w;
                dy_P_N  = y_N - y_P;   dx_P_N  = x_N - x_P;
                dy_N_Nw = y_Nw - y_N;  dx_N_Nw = x_Nw - x_N;
                dy_Nw_w = y_w - y_Nw;  dx_Nw_w = x_w - x_Nw;
                
                % Around etaw
                
                dy_W_P  = y_P - y_W;    dx_W_P  = x_P - x_W;
                dy_P_n  = y_n - y_P;    dx_P_n  = x_n - x_P;
                dy_n_nW = y_nW - y_n;   dx_n_nW = x_nW - x_n;
                dy_nW_W = y_W - y_nW;   dx_nW_W = x_W - x_nW;
                
                % Around P
                
                dy_n_nw = y_nw - y_n;  dx_n_nw = x_nw - x_n;
                dy_nw_w = y_w - y_nw;  dx_nw_w = x_w - x_nw;
                
                % At boundary
                
                dl_P_n = norm([dx_P_n; dy_P_n]);
                dl_w_P = norm([dx_w_P; dy_w_P]);
                
                % Areas
                
                S_etaw              = abs( x_n*y_P - y_n*x_P...
                    + x_P*y_W - y_P*x_W...
                    + x_W*y_nW - y_W*x_nW...
                    + x_nW*y_n - y_nW*x_n) / 2;
                S_nomega            = abs( x_N*y_P - y_N*x_P...
                    + x_P*y_w - y_P*x_w...
                    + x_w*y_Nw - y_w*x_Nw...
                    + x_Nw*y_N - y_Nw*x_N ) / 2;
                S_etaomega         = abs( x_n*y_P - y_n*x_P...
                    + x_P*y_w - y_P*x_w...
                    + x_w*y_nw - y_w*x_nw...
                    + x_nw*y_n - y_nw*x_n) / 2;
                
                
                %$$$$$$$$$$$$$$$$$$$$$$ Stecil $$$$$$$$$$$$$$$$$$$
                
                build_southEast
                
                %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                
                
                % West
                stecil(index(i, j-1)) = D_3;
                % North
                stecil(index(i-1, j)) = D_1;
                % NW
                stecil(index(i-1, j-1)) = D_4;
                
                
            end
            
            if strcmp(boundary.east, 'Neumann')
                % P
                stecil(index(i, j)) = D0;
                b = -(dl_P_n+dl_w_P)/S_etaomega*(-beta);
                
            elseif strcmp(boundary.east, 'Robin')
                
                % P
                stecil(index(i, j)) = D0 + (dl_P_n+dl_w_P)/S_etaomega*(-alpha/k);
                
                b = -(dl_P_n+dl_w_P)/S_etaomega*(-alpha/k*(-Tinf));
                
            elseif strcmp(boundary.east, 'Dirichlet')
                
                stecil(index(i, j)) = 1;
                b=TD.east;
            end
            
            
            
        elseif strcmp(boundary.south, 'Dirichlet')
            stecil(index(i, j)) = 1;
            b=TD.south;
        end
        
        
    case 'EastNorth'
        
        
        if strcmp(boundary.north, 'Neumann') || strcmp(boundary.north, 'Robin')
            
            if strcmp(boundary.east, 'Neumann') || strcmp(boundary.east, 'Robin')
                
                % Nomencature:
                %
                %  W=(i,j-1) -  -  w---- omega  ---P=(i,j)
                %                  |       |        |
                % sigmaW      sigmaw  sigmaomega   sigma
                %                  |       |        |
                %     sW - - - - -sw - -somega ---- s
                %                  |       |        |
                %                  |       |        |
                %                  |       |        |
                %    SW=(i+1,j-1)  Sw - Somega  -  S=(i+1,j)
                
                
                % Calculate some help values
                
                y_W  = Z(i  ,j-1);   x_W  = S(i  ,j-1);
                y_P  = Z(i  ,j  );   x_P  = S(i  ,j  );
                y_SW = Z(i+1,j-1);   x_SW = S(i+1,j-1);
                y_S  = Z(i+1,j  );   x_S  = S(i+1,j  );
                
                y_w  = (y_P  + y_W)/2;  x_w  = (x_P  + x_W)/2;
                y_Sw = (y_S  + y_SW)/2; x_Sw = (x_S  + x_SW)/2;
                y_sW = (y_W  + y_SW)/2;  x_sW = (x_W  + x_SW)/2;
                y_s  =  (y_P +  y_S)/2;  x_s  =  (x_P +  x_S)/2;
                y_sw = (y_w + y_Sw)/2;  x_sw = (x_w + x_Sw)/2;
                
                
                
                %{
                y_NE = Z(i-1,j+1);   x_NE = S(i-1,j+1);
                y_N  = Z(i-1,j  );   x_N  = S(i-1,j  );
                y_E  = Z(i  ,j+1);   x_E  = S(i  ,j+1);
                y_P  = Z(i  ,j  );   x_P  = S(i  ,j  );
                y_SE = Z(i+1,j+1);   x_SE = S(i+1,j+1);
                y_S  = Z(i+1,j  );   x_S  = S(i+1,j  );
                
                y_Ne = (y_N + y_NE)/2;  x_Ne = (x_N + x_NE)/2;
                y_e  = (y_P  + y_E)/2;  x_e  = (x_P  + x_E)/2;
                y_Se = (y_S  + y_SE)/2; x_Se = (x_S  + x_SE)/2;
                
                y_nE = (y_NE + y_E)/2;   x_nE = (x_NE + x_E)/2;
                y_sE = (y_E  + y_SE)/2;  x_sE = (x_E  + x_SE)/2;
                y_n  =  (y_N + y_P)/2;   x_n  =  (x_N + x_P)/2;
                y_s  =  (y_P +  y_S)/2;  x_s  =  (x_P +  x_S)/2;
                
                y_se = (y_e + y_Se)/2;  x_se = (x_e + x_Se)/2;
                y_ne = (y_e + y_Ne)/2;  x_ne = (x_e + x_Ne)/2;
                %}
                
                % Around somega
                
                dy_Sw_S = y_S - y_Sw;  dx_Sw_S = x_S - x_Sw;
                dy_S_P  = y_P - y_S;   dx_S_P  = x_P - x_S;
                dy_P_w  = y_w - y_P;   dx_P_w  = x_w - x_P;
                dy_w_Sw = y_Sw - y_w;  dx_w_Sw = x_Sw - x_w;
                
                % Around sigmaw
                
                dy_sW_s = y_s - y_sW;  dx_sW_s = x_s - x_sW;
                dy_s_P  = y_P - y_s;   dx_s_P  = x_P - x_s;
                dy_P_W  = y_W - y_P;   dx_P_W  = x_W - x_P;
                dy_W_sW = y_sW - y_W;  dx_W_sW = x_sW - x_W;
                
                % Around P
                
                dy_sw_s  = y_s - y_sw;  dx_sw_s  = x_s - x_sw;
                dy_w_sw  = y_sw - y_w;  dx_w_sw  = x_sw - x_w;
                
                % At boundary
                
                dl_s_P = norm([dx_s_P; dy_s_P]);
                dl_P_w = norm([dx_P_w; dy_P_w]);
                
                % Areas
                
                S_sigmaomega          = abs( x_P*y_s - y_P*x_s...
                    + x_s*y_sw - y_s*x_sw...
                    + x_sw*y_w - y_sw*x_w...
                    + x_w*y_P - y_w*x_P) / 2;
                S_somega                 = abs( x_P*y_S - y_P*x_S...
                    + x_S*y_Sw - y_S*x_Sw...
                    + x_Sw*y_w - y_Sw*x_w...
                    + x_w*y_P - y_w*x_P) / 2;
                S_sigmaw                 = abs( x_P*y_s - y_P*x_s...
                    + x_s*y_sW - y_s*x_sW...
                    + x_sW*y_W - y_sW*x_W...
                    + x_W*y_P - y_W*x_P) / 2;
                
                
                %$$$$$$$$$$$$$$$$$$$$$$ Stecil $$$$$$$$$$$$$$$$$$$
                
                build_northEast
                
                %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                
                
                % West
                stecil(index(i, j-1)) = D_3;
                % South
                stecil(index(i+1, j)) = D1;
                % SW
                stecil(index(i+1, j-1)) = D_2;
                
                
            end
            
            if strcmp(boundary.east, 'Neumann')
                % P
                stecil(index(i, j)) = D0;
                b = -(dl_s_P+dl_P_w)/S_sigmaomega*(-beta);
                
            elseif strcmp(boundary.east, 'Robin')
                
                % P
                stecil(index(i, j)) = D0 + (dl_s_P+dl_P_w)/S_sigmaomega*(-alpha/k);
                
                b = -(dl_s_P+dl_P_w)/S_sigmaomega*(-alpha/k*(-Tinf));
                
            elseif strcmp(boundary.east, 'Dirichlet')
                
                stecil(index(i, j)) = 1;
                b=TD.east;
            end
            
            
            
            %        stecil(index(i,j)) = 1;
            %        if i == 1 && ~(strcmp(boundary.north, 'Dirichlet') && strcmp(boundary.east, 'Dirichlet'))
            %            stecil(index(i, j-1)) = -0.5;
            %            stecil(index(i+1, j)) = -0.5;
            %        else
            %            stecil(index(i-1, j)) = -1;
            %        end
            
            
        elseif strcmp(boundary.north, 'Dirichlet')
            stecil(index(i, j)) = 1;
            b=TD.north;
        end
        
        
    case 'WestSouth'
        
        if strcmp(boundary.south, 'Neumann')
            
            if strcmp(boundary.west, 'Neumann') || strcmp(boundary.west, 'Robin')
                
                
                
                % Calculate some help values
                
                y_NE = Z(i-1,j+1);   x_NE = S(i-1,j+1);
                y_N  = Z(i-1,j  );   x_N  = S(i-1,j  );
                y_E  = Z(i  ,j+1);   x_E  = S(i  ,j+1);
                y_P  = Z(i  ,j  );   x_P  = S(i  ,j  );
                
                y_Ne = (y_N + y_NE)/2;  x_Ne = (x_N + x_NE)/2;
                y_e  = (y_P  + y_E)/2;  x_e  = (x_P  + x_E)/2;
                
                y_nE = (y_NE + y_E)/2;   x_nE = (x_NE + x_E)/2;
                y_n  =  (y_N + y_P)/2;   x_n  =  (x_N + x_P)/2;
                
                y_ne = (y_e + y_Ne)/2;  x_ne = (x_e + x_Ne)/2;
                
                % Around nepsilon
                
                dy_P_e   = y_e - y_P;    dx_P_e   = x_e - x_P;
                dy_e_Ne  = y_Ne - y_e;   dx_e_Ne  = x_Ne - x_e;
                dy_Ne_N  = y_N - y_Ne;   dx_Ne_N  = x_N - x_Ne;
                dy_N_P   = y_P - y_N;    dx_N_P   = x_P - x_N;
                
                
                % Around etae
                
                dy_P_E  = y_E - y_P;    dx_P_E  = x_E - x_P;
                dy_E_nE = y_nE - y_E;   dx_E_nE = x_nE - x_E;
                dy_nE_n = y_n - y_nE;   dx_nE_n = x_n - x_nE;
                dy_n_P  = y_P - y_n;    dx_n_P  = x_P - x_n;
                
                
                % Around P
                
                dy_e_ne = y_ne - y_e;  dx_e_ne = x_ne - x_e;
                dy_ne_n = y_n - y_ne;  dx_ne_n = x_n - x_ne;
                
                % At boundary
                
                dl_n_P = norm([dx_n_P; dy_n_P]);
                dl_P_e = norm([dx_P_e; dy_P_e]);
                
                % Areas
                
                S_etae                 = abs( x_nE*y_E - y_nE*x_E...
                    + x_E*y_P - y_E*x_P...
                    + x_P*y_n - y_P*x_n...
                    + x_n*y_nE - y_n*x_nE) / 2;
                
                
                S_nepsilon                 = abs( x_P*y_N - y_P*x_N...
                    + x_N*y_Ne - y_N*x_Ne...
                    + x_Ne*y_e - y_Ne*x_e...
                    + x_e*y_P - y_e*x_P ) / 2;
                
                
                S_etaepsilon         = abs( x_P*y_n - y_P*x_n...
                    + x_n*y_ne - y_n*x_ne...
                    + x_ne*y_e - y_ne*x_e...
                    + x_e*y_P - y_e*x_P) / 2;
                
                
                
                
                %$$$$$$$$$$$$$$$$$$$$$$ Stecil $$$$$$$$$$$$$$$$$$$
                
                build_southWest
                
                %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                
                % East
                stecil(index(i, j+1)) = D3;
                % North
                stecil(index(i-1, j)) = D_1;
                % NE
                stecil(index(i-1, j+1)) = D2;
                
                
            end
            
            %%% NOTE: There is one face missing for the stencil in P (dl_P_e)
            
            if strcmp(boundary.west, 'Neumann')
                % P
                stecil(index(i, j)) = D0;
                b = -(dl_n_P+dl_P_e)/S_etaepsilon*(-beta);
                
            elseif strcmp(boundary.west, 'Robin')
                
                % P
                stecil(index(i, j)) = D0 + (dl_n_P+dl_P_e)/S_etaepsilon*(-alpha/k);
                
                b = -(dl_n_P+dl_P_e)/S_etaepsilon*(-alpha/k*(-Tinf));
                
            elseif strcmp(boundary.west, 'Dirichlet')
                
                stecil(index(i, j)) = 1;
                b=TD.west;
            end
            
            
            
        elseif strcmp(boundary.south, 'Dirichlet')
            stecil(index(i, j)) = 1;
            b=TD.south;
        end
        
        
    case 'WestNorth'
        
        
        if strcmp(boundary.north, 'Neumann') || strcmp(boundary.north, 'Robin')
            
            if strcmp(boundary.west, 'Neumann') || strcmp(boundary.west, 'Robin')
                

                
                
                % Calculate some help values
                
                y_E  = Z(i  ,j+1);   x_E  = S(i  ,j+1);
                y_P  = Z(i  ,j  );   x_P  = S(i  ,j  );
                y_SE = Z(i+1,j+1);   x_SE = S(i+1,j+1);
                y_S  = Z(i+1,j  );   x_S  = S(i+1,j  );
                
                y_e  = (y_P  + y_E)/2;  x_e  = (x_P  + x_E)/2;
                y_Se = (y_S  + y_SE)/2; x_Se = (x_S  + x_SE)/2;
                y_sE = (y_E  + y_SE)/2;  x_sE = (x_E  + x_SE)/2;
                y_s  =  (y_P +  y_S)/2;  x_s  =  (x_P +  x_S)/2;
                y_se = (y_e + y_Se)/2;  x_se = (x_e + x_Se)/2;
                
                % Around sepsilon
                
                dy_S_Se  = y_Se - y_S;   dx_S_Se = x_Se - x_S;
                dy_Se_e  = y_e - y_Se;   dx_Se_e = x_e - x_Se;
                dy_e_P   = y_P - y_e;    dx_e_P  = x_P - x_e;
                dy_P_S   = y_S - y_P;    dx_P_S  = x_S - x_P;
                
                
                % Around sigmae
                
                dy_s_sE  = y_sE - y_s ; dx_s_sE  = x_sE - x_s ;
                dy_sE_E  = y_E - y_sE ; dx_sE_E  = x_E - x_sE ;
                dy_E_P   = y_P - y_E ;  dx_E_P   = x_P - x_E ;
                dy_P_s   = y_s - y_P ;  dx_P_s   = x_s - x_P ;
                
                
                % Around P
                
                dy_s_se  = y_se - y_s;  dx_s_se  = x_se - x_s;
                dy_se_e  = y_e - y_se;  dx_se_e  = x_e - x_se;
                
                % At boundary
                
                dl_P_s = norm([dx_P_s; dy_P_s]);
                dl_e_P = norm([dx_e_P; dy_e_P]);
                
                % Areas
                
                
                S_sigmaepsilon  = abs( x_s*y_P - y_s*x_P...
                    + x_P*y_e - y_P*x_e...
                    + x_e*y_se - y_e*x_se...
                    + x_se*y_s - y_se*x_s) / 2;
                
                
                S_sepsilon                 = abs( x_S*y_P - y_S*x_P...
                    + x_P*y_e - y_P*x_e...
                    + x_e*y_Se - y_e*x_Se...
                    + x_Se*y_S - y_Se*x_S) / 2;
                
                
                S_sigmae                 = abs( x_E*y_sE - y_E*x_sE...
                    + x_sE*y_s - y_sE*x_s...
                    + x_s*y_P - y_s*x_P...
                    + x_P*y_E - y_P*x_E) / 2;
                
                
                %$$$$$$$$$$$$$$$$$$$$$$ Stecil $$$$$$$$$$$$$$$$$$$
                
                build_northWest
                
                %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                
                % East
                stecil(index(i, j+1)) = D3;
                % South
                stecil(index(i+1, j)) = D1;
                % SE
                stecil(index(i+1, j+1)) = D4;
                
                
            end
            
            if strcmp(boundary.west, 'Neumann')
                % P
                stecil(index(i, j)) = D0;
                b = -(dl_P_s+dl_e_P)/S_sigmaepsilon*(-beta);
                
            elseif strcmp(boundary.west, 'Robin')
                
                % P
                stecil(index(i, j)) = D0 + (dl_P_s+dl_e_P)/S_sigmaepsilon*(-alpha/k);
                
                b = -(dl_P_s+dl_e_P)/S_sigmaepsilon*(-alpha/k*(-Tinf));
                
            elseif strcmp(boundary.west, 'Dirichlet')
                
                stecil(index(i, j)) = 1;
                b=TD.west;
                
            end
            
            
            
            %        stecil(index(i,j)) = 1;
            %        if i == 1 && ~(strcmp(boundary.north, 'Dirichlet') && strcmp(boundary.east, 'Dirichlet'))
            %            stecil(index(i, j-1)) = -0.5;
            %            stecil(index(i+1, j)) = -0.5;
            %        else
            %            stecil(index(i-1, j)) = -1;
            %        end
            
            
        elseif strcmp(boundary.north, 'Dirichlet')
            stecil(index(i, j)) = 1;
            b=TD.north;
        end
        
        
end



