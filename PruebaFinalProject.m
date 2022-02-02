% Numerical code to solve the 2D steady heat equation by Finite Differences
% Second order accurate in Space
% Dirichlet West, Neumann South, Robin North & East


%% Cleaning %%
close all
clear
clc
format shortG
tic
addpath(genpath('GenerateStencil'))

PruebainitDat
%------------------------------------------------------
%% 2D MESH (Matrix S and Matrix Y)

dX= l/(dimX-1);
% dY= h/(dimY-1);
dZ= f/(dimZ-1);

s= linspace(0,l,dimX);
% y= linspace(h,0,dimY);
z= linspace(f,0,dimZ);
zf= z(1,1:z2);
zb= z(1,z2+1:dimZ);

S= repmat(s,[dimZ,1]);
% Y= repmat(y',[1,dimX]);
Z= repmat(z',[1,dimX]);
Zf= repmat(zf',[1,dimX]);
Zb= repmat(zb',[1,dimX]);

%% Temporal dimensions

tend= 1;
dt = 0.01*1/(2*D(1))*(dX.^2.*dZ.^2)/(dX.^2+dZ.^2);
% dt for weighted average
% theta= 0.5
% dt = 1/(2*D(1)*(1-2*theta*(theta<0.5)))*(dX^2.*dZ^2)/(dX^2+dZ^2);
t = 0:dt:tend;
tn= numel(t)

%% Defining the equation (A·T=b) %%
%For the board
A1= zeros(dimXz1,dimXz1);
B1= zeros(dimXz1,1);
Tb0= zeros(dimXz1,1);      %Temperature array

%For the fluid
A2= zeros(dimXz2,dimXz2);
B2= zeros(dimXz2,1);
Tf0= zeros(dimXz2,1);      %Temperature array


%% Thermal conductivity  matrix K (nodal) %%

K= ones(dimZ,dimX);
K(1:z2,:)= K(1:z2,:)*Kval(2);
K(z2+1:dimZ,:)= K(z2+1:dimZ,:)*Kval(1);
Kb= ones(z1,dimX)*Kval(1);
Kf= ones(z2,dimX)*Kval(2);
% if strcmp(heat_conduc,'non homogeneous')
%     flag= (S>=xk(1)).*(S<=xk(2)).*(Y>=yk(1)).*(Y<=yk(2));
%     K(flag==1)= KnH;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Finite Volumes
% Index maps the node position to the correct linear equation
index = @(ii, jj) ii + (jj-1) * z1;
% set up the system matrix A
Kval= Kb(1,1);
for i = 1:z1
    for j = 1:dimX
        % Fill the system matrix and the right-hand side for node (i,j)
        [A1(index(i,j), :),B1(index(i,j))] = ...
            PruebastampFP(i, j, S(z2+1:end,:), Zb, TD, alpha, Tinf, boundary, Kval, beta);
    end
end

% BC and temperatures for fluid
boundary.north = 'Dirichlet';   TD.north = 298;
boundary.south = 'Dirichlet';   TD.south= 298;
boundary.west = 'Dirichlet';    TD.west= 298;
boundary.east = 'Robin';    TD.east= 298;

% set up the system matrix A2
Kval= Kf(1,1);
for i = 1:z2
    for j = 1:dimX
        % Fill the system matrix and the right-hand side for node (i,j)
        [A2(index(i,j), :),B2(index(i,j))] = ...
            PruebastampFP(i, j, S(z2+1:end,:), Zf, TD, alpha, Tinf, boundary, Kval, beta);
    end
end

% CONVECTION
% Inner nodes for fluid A2%
A3=zeros(size(A2));
for i= 2:dimX-1
    for j= 2:z2-1
        ms= z2*(i-1) + j;  %counter used for the loop

        A2(ms,ms)= A2(ms,ms) - v_x(j)*(1/dX);  %Position component
        A2(ms,ms-z1)= A2(ms,ms-z1) + v_x(j)*1/dX;       %West component
    end
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% Defining the Source %
% % Fluid
% WF= zeros(size(B2));
% if strcmp(source,'yes')
%     flag= ( ((S(1:z2,:)-c3).^2 + (Zf-d3).^2) <= R3.^2 );
%     WF(flag==1)= w3;
%     loc= find(WF);
%     for k = 1:length(loc)
%      A2(loc(k),:)= 0;
%      A2(loc(k),loc(k))= 1;
%      B2(loc(k))= WF(loc(k));
%      Tnf = Tnf(:);
%      Tnf(loc(k))= w3;
%      Tnf = reshape(Tnf,z2,dimX);
%     end
% %     flag= (Z>=z1*dZ).*(Z<=f);
% %     W(flag==1)= 0;
% end
%  WFindex= reshape(WF,z2,dimX);  % for the plotting

% Board
WB= zeros(size(B1));

if strcmp(source,'yes')
    flag= ( ((S(z2+1:dimZ,:)-c1).^2 + (Zb-d1).^2) <= R1.^2 );
    WB(flag==1)= w1;
    loc= find(WB);
    for k = 1:length(loc)
     A1(loc(k),:)= 0;
     A1(loc(k),loc(k))= 1;
     B1(loc(k))= WB(loc(k));
     Tnb = Tnb(:);
     Tnb(loc(k))= w1;
     Tnb = reshape(Tnb,z1,dimX);
    end
    flag= (Z>=z1*dZ).*(Z<=f);
    W(flag==1)= 0;
end

 WBindex= reshape(WB,z1,dimX);  % for the plotting
 
% % %  %% Source 2
% % %  
% % % WB= zeros(size(B1));
% % % if strcmp(source,'yes')
% % %     flag= (((S(z2+1:dimZ,:)-c2).^2 + (Zb-d2).^2) <= R2.^2 );
% % %     WB(flag==1)= w2;
% % %     loc= find(WB);
% % %     for k = 1:length(loc)
% % %      A1(loc(k),:)= 0;
% % %      A1(loc(k),loc(k))= 1;
% % %      B1(loc(k))= WB(loc(k));
% % %      Tnb = Tnb(:);
% % %      Tnb(loc(k))= w2;
% % %      Tnb = reshape(Tnb,z1,dimX);
% % %     end
% % %     flag= (Z>=z1*dZ).*(Z<=f);
% % %     W(flag==1)= 0;
% % % end
% % % 
% % % WB= zeros(size(B1));
% % % 
% % % third source
% % % for c=5:15
% % %     WB(index(z2+4,c))=w2;
% % % end
% % % loc= find(WB);
% % % for k = 1:length(loc)
% % %  A1(loc(k),:)= 0;
% % %  A1(loc(k),loc(k))= 1;
% % %  B1(loc(k))= WB(loc(k));
% % %  Tnb = Tnb(:);
% % %  Tnb(loc(k))= w2;
% % %  Tnb = reshape(Tnb,z1,dimX);
% % % end
% % % 
% % %  
 WBindex= WBindex+reshape(WB,z1,dimX);  % for the plotting
 B1= reshape(B1,z1,dimX);
 B2= reshape(B2,z2,dimX);
%-----------------------------------------------------------------------
A1= sparse(A1);
A2= sparse(A2);
B1= sparse(B1);
B2= sparse(B2);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% Temporal evolution

Tt= zeros(dimZ,dimX,tn);

Tnf(:,end) = TD.east*strcmp(boundary.east,'Dirichlet')+...
    (1-strcmp(boundary.east,'Dirichlet'))*Tnf(:,end);
Tnf(end,:) = TD.south*strcmp(boundary.south,'Dirichlet')+...
    (1-strcmp(boundary.south,'Dirichlet'))*Tnf(end,:);
Tnf(:,1) = TD.west*strcmp(boundary.west,'Dirichlet')+...
    (1-strcmp(boundary.west,'Dirichlet'))*Tnf(:,1);
Tnf(1,:) = TD.north*strcmp(boundary.north,'Dirichlet')+...
    (1-strcmp(boundary.north,'Dirichlet'))*Tnf(1,:);
%% Runge-Kutta 4

tic
for i=1:tn-1
    [tn i]
    
    Tp_b = Tnb(:) + dt/2*A1*Tnb(:) - B1(:)*dt/2;
    Tpp_b = Tnb(:) + dt/2*A1*Tp_b(:) - B1(:)*dt/2;
    Tppp_b = Tnb(:) + dt*A1*Tpp_b(:) - B1(:)*dt;
    Tb(:,:,i+1) = reshape(Tnb(:)+ 1/6*dt*(A1*Tnb(:) +2*A1*Tp_b(:)+...
           2*A1*Tpp_b(:) +A1*Tppp_b(:))- B1(:)*dt, z1,dimX);
    Tnb = Tb(:,:,i+1);  % For the next loop
    
    Tnf(end,:) = Tb(1,:,i+1);  %%solo valido para dimensiones iguales
    B2(end,:) = Tb(1,:,i+1);
    Tp_f = Tnf(:) + dt/2*A2*Tnf(:) - B2(:)*dt/2;
    Tpp_f = Tnf(:) + dt/2*A2*Tp_f(:) - B2(:)*dt/2;
    Tppp_f = Tnf(:) + dt*A2*Tpp_f(:) - B2(:)*dt;
    Tf(:,:,i+1) = reshape(Tnf(:)+1/6*dt*(A2*Tnf(:)+2*A2*Tp_f(:)...
           +2*A2*Tpp_f(:)+A2*Tppp_f(:))-B2(:)*dt,z2,dimX);
    Tnf = Tf(:,:,i+1);
    
    Tt(1:z2,:,i)= Tf(:,:,i);  % Total temperature matrix (plots)
    Tt(z2+1:dimZ,:,i)= Tb(:,:,i);
    
end
elapsed_t= toc
tplot= (1:tn)*dt;  % For the plots
tic
%     Tt(1:z2,:,:)= Tf(:,:,:);  % Total temperature matrix (plots)
%     Tt(z2+1:dimZ,:,:)= Tb(:,:,:);
toc
%% Ploting Results %%
% 
% %% Total
% figure
% contour(S,Z,Tt(:,:,tn-1),150)
% title('Temperature distribution (contour)')
% % title(['Temperature distribution (contour), \alpha=',num2str(alpha)])
%     xlabel('x')
%     ylabel('z')
%     colorbar
%     xlim([0, l]);
%     ylim([0, f]);
%     print(gcf,'LEDnoconv','-dpng')
 
%     %%
% figure
% surf(S,Z,Tt(:,:,tn-1));view(2);
% title('Total T distribution (surface)')
%     xlabel('x')
%     ylabel('z')
%     colorbar
%     caxis([250 350]);
%     xlim([0, l]);
%     ylim([0, f]); 
%     print(gcf,'LEDnoconv','-dpng')
%% GIFs
figure
  filename = '2LEDconv.gif';
for k=1
    pcolor(S,Z,Tt(:,:,k))
    colorbar
%    title(['dt=',num2str(dt)])
    drawnow
    frame=getframe(gcf);
    im=frame2im(frame);
    [imind,cm]=rgb2ind(im,8);
    imwrite(imind,cm,filename,'gif','Loopcount',inf)
end
for k=2:10000:tn
    pcolor(S,Z,Tt(:,:,k))
    colorbar
    caxis([298 350]);
    colorbar    
    hold on
    yline(1, 'r','LineWidth',1);
    hold on
%    title(['Board dt= ',num2str(dt),', t= ',num2str(tplot(k))])
    title(['Distribution of T, t= ',num2str(tplot(k)),'s, \alpha=',num2str(alpha)])
    drawnow
    frame=getframe(gcf);
    im=frame2im(frame);
    [imind,cm]=rgb2ind(im,8);
    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0)
    k;
end
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Filled contour
% Tt= zeros(dimZ,dimX,tn);
figure
filename='1LEDfilled.gif';
for k=1
    Tt(1:z2,:,k)= Tf(:,:,k);
    Tt(z2+1:dimZ,:,k)= Tb(:,:,k);
    contourf(S,Z,Tt(:,:,k))
    hold on 
    title(['dt=',num2str(dt)])
    drawnow
    frame=getframe(gcf);
    im=frame2im(frame);
    [imind,cm]=rgb2ind(im,8);
    imwrite(imind,cm,filename,'gif','Loopcount',inf)
end
for k=2:1000:tn-1
    Tt(1:z2,:,k)= Tf(:,:,k);
    Tt(z2+1:dimZ,:,k)= Tb(:,:,k);
    hold on
    yline(1, 'r','LineWidth',1);
    hold on
    contourf(S,Z,Tt(:,:,k));
    caxis([298 350]);
    colorbar
    title(['Board dt= ',num2str(dt),', t= ',num2str(tplot(k))])
    title(['Isotherms, t= ',num2str(tplot(k)),'s, \alpha=',num2str(alpha)])
    drawnow
    frame=getframe(gcf);
    im=frame2im(frame);
    [imind,cm]=rgb2ind(im,8);
    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0)
    k;
end
%end

% %% Fluid Contour plot 
% figure
% contour(S(1:z2,:),Zf,Tf(:,:,tn-1),150)
% % title('Temperature distribution (contour)')
% title(['Temperature distribution (contour), \alpha=',num2str(alpha)])
%     xlabel('x')
%     ylabel('z')
%     colorbar
%     xlim([0, l]);
%     ylim([z1, f]);
%     print(gcf,'contourSpecial','-dpng')
    
 %% Board RK
% tic 
% figure
%   filename = 'fvm_mainf.gif';
% for k=1
%     pcolor(S(z2+1:dimZ,:),Zb,Tb(:,:,k))
%     colorbar
%     title(['dt=',num2str(dt)])
%     drawnow
%     frame=getframe(gcf);
%     im=frame2im(frame);
%     [imind,cm]=rgb2ind(im,8);
%     imwrite(imind,cm,filename,'gif','Loopcount',inf)
% end
% for k=2:1000:tn
%     pcolor(S(z2+1:dimZ,:),Zb,Tb(:,:,k))
%     colorbar
%     caxis([298 398]); % Change depending on T values selected
%     title(['Board dt= ',num2str(dt),', t= ',num2str(tplot(k))])
%     drawnow
%     frame=getframe(gcf);
%     im=frame2im(frame);
%     [imind,cm]=rgb2ind(im,8);
%     imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0)
%     k;
% end
% toc
% 
% %% Fluid RK
% tic 
% figure
%   filename = 'fluidSpecial.gif';
% for k=1
%     pcolor(S(1:z2,:),Zf,Tf(:,:,k))
%     colorbar
%     title(['dt=',num2str(dt)])
%     drawnow
%     frame=getframe(gcf);
%     im=frame2im(frame);
%     [imind,cm]=rgb2ind(im,8);
%     imwrite(imind,cm,filename,'gif','Loopcount',inf)
% end
% for k=2:1000:tn
%     pcolor(S(1:z2,:),Zf,Tf(:,:,k))
%     colorbar
%     caxis([298 398]);
%     colorbar
% %     title(['Board dt= ',num2str(dt),', t= ',num2str(tplot(k))])
%     title(['Fluid dt= ',num2str(dt),', t= ',num2str(tplot(k)),'s, \alpha=',num2str(alpha)])
%     drawnow
%     frame=getframe(gcf);
%     im=frame2im(frame);
%     [imind,cm]=rgb2ind(im,8);
%     imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0)
%     k;
% end
% toc   
%% 
% %%% K and Source %%% 
% % Mesh plot: Thermal conductivity
% figure
% % subplot(121)
% mesh(S,Z,K); view(2);
% title('Thermal conductivity distribution (mesh)')
%     xlabel('x')
%     ylabel('z')
%     colorbar
%     xlim([0, l]);
%     ylim([0, f]);
%     print(gcf,'Kdistr','-dpng')
% figure    
% % Mesh plot: pointwise source  
% % subplot(122)
% mesh(S(z2+1:dimZ,:),Zb,WBindex); view(2);
% axis equal   %to notice the circle
% title('Pointwise source distribution (mesh)')
%     xlabel('x')
%     ylabel('z')
%     xlim([0, l]);
%     ylim([0, f]);
%     colorbar
%      print(gcf,'Source','-dpng')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
  

% Tt= zeros(dimZ,dimX);
% Tt(1:z2,:)= Tf(:,:,end);
% Tt(z2+1:dimZ,:)= Tb(:,:,end);


% 
% % Contour plot
% subplot(122)
% contour(S,Z,Tt,200)
% title('Total T distribution (contour)')
%     xlabel('x')
%     ylabel('z')
%     colorbar
%     xlim([0, l]);
%     ylim([0, f]);
% % print(gcf,'2DHeatEq','-dpng')
% 

% %% Board %%%
% % Surface plot
% figure (3)
% subplot(221)
% surf(S(z2+1:dimZ,:),Zb,Tb0);view(2);
% title('Board: T distribution (surface)')
%     xlabel('x')
%     ylabel('z')
%     colorbar
%     xlim([0, l]);
%     ylim([0, (z1-1)*dZ]);
% 
% % Contour plot
% subplot(222)
% contour(S(z2+1:dimZ,:),Zb,Tb0,150)
% title('Temperature distribution (contour)')
%     xlabel('x')
%     ylabel('z')
%     colorbar
%     xlim([0, l]);
%     ylim([0, (z1-1)*dZ]);
%     
% 
% % print(gcf,'2DHeatEq','-dpng')
% 
% Tt= zeros(dimZ,dimX,tn);
% figure
% %% surface
% filename='surface.gif',
% for k=1
%     Tt(1:z2,:,k)= Tf(:,:,k);
%     Tt(z2+1:dimZ,:,k)= Tb(:,:,k);
%     pcolor(S,Z,Tt(:,:,k))
%     hold on 
%     caxis([0 120]);
%     title(['dt=',num2str(dt)])
%     drawnow
%     frame=getframe(gcf);
%     im=frame2im(frame);
%     [imind,cm]=rgb2ind(im,8);
%     imwrite(imind,cm,filename,'gif','Loopcount',inf)
% end
% for k=2:10000:tn
%     Tt(1:z2,:,k)= Tf(:,:,k);
%     Tt(z2+1:dimZ,:,k)= Tb(:,:,k);
%     pcolor(S,Z,Tt(:,:,k))
%     hold on 
%     yline(1, 'LineWidth',1);
%     caxis([10 100]);
%     colorbar
%     title(['Board dt= ',num2str(dt),', t= ',num2str(tplot(k)),'s'])
%     drawnow
%     frame=getframe(gcf);
%     im=frame2im(frame);
%     [imind,cm]=rgb2ind(im,8);
%     imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0)
%     k;
% end

%% Fluid %%%
% % Surface plot
% subplot(223)
% surf(S(1:z2,:),Zf,Tf0);view(2);
% title('Fluid: T distribution (surface)')
%     xlabel('x')
%     ylabel('z')
%     colorbar
%     xlim([0, l]);
%     ylim([z1*dZ, f]);
% 
% % Contour plot
% subplot(224)
% contour(S(1:z2,:),Zf,Tf(:,:,tn-1),150)
% title(['Temperature distribution (contour), \alpha=',num2str(alpha)])
%     xlabel('x')
%     ylabel('z')
%     colorbar
%     xlim([0, l]);
%     ylim([z1*dZ, f]);
%     print(gcf,'contourRobina1000','-dpng')    
% % % legend({'y = sin(x)','y = cos(x)'},'Location','southwest')
% % % title(gcf,'fluid')??
% % % k = sin(pi/2); title(['sin(\pi/2) = ' num2str(k)])

