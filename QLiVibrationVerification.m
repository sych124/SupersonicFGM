%................................................................
 
% 1st Paper Project in 2015 
% Vibaratino analysis of functionally graded plates in thermal effects
% Written by Kang

% Concept

% 1. Nonlinear/uniform temperature distribution
% 2. Voigt-rule 
% 3. P-FGM
% 4. Neutral plane
% 5. von-Karman nonlinear term
% 6. Based on first-order shear deformation plate theory


% clear memory
clear all;colordef white;
clc
format short

kt=[2];
T_m=300;
dT=[500];
for tk=1:size(dT,2);
% materials
poisson = 0.30; kapa=5/6; 
thickness=1*1;h=thickness;
k_index=kt(1);
a=10; % a=Lx
b=10; % b=Ly
T_c=T_m+dT(tk);

% Property of FGplate
n=100;
aLayer=h/n;

% For z0 of neutral plan
zA=0;
zB=0;
for i=-h/2:aLayer:h/2  % for z0
    
    % Uniform temperature
    Tzz=T_c;
%     % Nonlinear temperature(polynominal)
%     if i<-h/2+aLayer;
%     [E_c,E_m,a_c,a_m,Rho_c,Rho_m,K_c1,K_m1]=MaterialPropertyWithTemperature(T_m);
%     K_m=K_m1;
%     K_c=K_c1;
%     end
%     Kcm=K_c-K_m;
%     P=i/h+0.5;
%     C=1-Kcm/((k_index+1)*K_m)+Kcm^2/((2*k_index+1)*K_m^2)-Kcm^3/((3*k_index+1)*K_m^3)+Kcm^4/((4*k_index+1)*K_m^4)...
%         -Kcm^5/((5*k_index+1)*K_m^5);
%     Rz=(P-Kcm/((k_index+1)*K_m)*P^(k_index+1)+Kcm^2/((2*k_index+1)*K_m^2)*P^(2*k_index+1)...
%         -Kcm^3/((3*k_index+1)*K_m^3)*P^(3*k_index+1)+Kcm^4/((4*k_index+1)*K_m^4)*P^(4*k_index+1)...
%         -Kcm^5/((5*k_index+1)*K_m^5)*P^(5*k_index+1))/C;
%     Tzz=T_m+(T_c-T_m)*Rz;

    [E_c,E_m,a_c,a_m,Rho_c,Rho_m,K_c,K_m]=MaterialPropertyWithTemperature(Tzz);

    Ez0=E_m+(E_c-E_m)*(i/h+0.5)^k_index;  % E(z) of power law
    Rhoz=Rho_m+(Rho_c-Rho_m)*(i/h+0.5)^k_index;  % Rho(z) of power law
    E=Ez0/(1-poisson^2)*[1 poisson 0;poisson 1 0; 0 0 (1-poisson)/2];

    zA=zA+Ez0*i;
    zB=zB+Ez0;
        
end
z0=zA/zB;
%z0=0;

% For ABD matrix
A= zeros(3,3);
B=zeros(3,3);
D=zeros(3,3);
Ss=zeros(2,2);
I0=0;
I1=0;
I2=0;
N_T=zeros(3,1);

% To get the ABD matrics
zw=zeros(size(-h/2:aLayer:h/2,2),1);
zw(:,1)=-h/2:aLayer:h/2;

for i=1:size(zw,1)   % for ABD matrix
    
    % Uniform temperature
    Tz(i,tk)=T_c;       
%     % Nonlinear temperature(polynominal)
%     if i<2;
%     [E_c,E_m,a_c,a_m,Rho_c,Rho_m,K_c1,K_m1]=MaterialPropertyWithTemperature(T_m);
%     K_m=K_m1;
%     K_c=K_c1;
%     end
%     Kcm=K_c-K_m;
%     P=zw(i)/h+0.5;
%     C=1-Kcm/((k_index+1)*K_m)+Kcm^2/((2*k_index+1)*K_m^2)-Kcm^3/((3*k_index+1)*K_m^3)+Kcm^4/((4*k_index+1)*K_m^4)...
%         -Kcm^5/((5*k_index+1)*K_m^5);
%     Rz=(P-Kcm/((k_index+1)*K_m)*P^(k_index+1)+Kcm^2/((2*k_index+1)*K_m^2)*P^(2*k_index+1)...
%         -Kcm^3/((3*k_index+1)*K_m^3)*P^(3*k_index+1)+Kcm^4/((4*k_index+1)*K_m^4)*P^(4*k_index+1)...
%         -Kcm^5/((5*k_index+1)*K_m^5)*P^(5*k_index+1))/C;
%     Tz(i,tk)=T_m+(T_c-T_m)*Rz;
    
    % Material properties in thermal effect each layer
    [E_c,E_m,a_c,a_m,Rho_c,Rho_m,K_c,K_m,v_c,v_m]=MaterialPropertyWithTemperature(Tz(i,tk));
    
    Ez(i,tk)=E_m+(E_c-E_m)*(zw(i)/h+0.5)^k_index;  % E(z) of power law
    Rhoz(i,tk)=Rho_m+(Rho_c-Rho_m)*(zw(i)/h+0.5)^k_index;  % Rho(z) of power law
    az(i,tk)=a_m+(a_c-a_m)*(zw(i)/h+0.5)^k_index;  % Rho(z) of power law
    E=Ez(i,tk)/(1-poisson^2)*[1 poisson 0;poisson 1 0; 0 0 (1-poisson)/2];
    Ni=zw(i)-z0;
    A=A+E*aLayer;
    B=B+E*Ni*aLayer;
    D=D+E*Ni^2*aLayer;
    Ss=Ss+1/2/(1+poisson)*Ez(i,tk)*[1 0;0 1]*aLayer*kapa;
    I0=I0+Rhoz(i,tk)*aLayer;
    I1=I1+Rhoz(i,tk)*Ni*aLayer;
    I2=I2+Rhoz(i,tk)*Ni^(2)*aLayer;
    N_T=N_T+E*[az(i,tk);az(i,tk);az(i,tk)]*aLayer*(Tz(i,tk)-T_m);
    %N_T(i,tk)=N_T(i,tk)+E*az(i,tk)*aLayer*(Tz(i,tk)-T_m);
    %M_T=M_T+E*[az(i,tk);az(i,tk);az(i,tk)]*aLayer*Ni*(Tz(i,tk)-T_m);
  
end

%Mesh generation
% Generation of coordinates(9 node)
Lx=a;Ly=b;
numberElementsX=4;
numberElementsY=4;
numberElements=numberElementsX*numberElementsY;
nodeCoordinatesX = linspace(0,Lx,2*numberElementsX+1)';
nodeCoordinatesY = linspace(0,Ly,2*numberElementsY+1)';
numberNodesX=size(nodeCoordinatesX,1);
numberNodesY=size(nodeCoordinatesY,1);
for i=1:numberNodesY
      for j=1:numberNodesX
        nodeCoordinates((i-1)*numberNodesX+j,1)=nodeCoordinatesX(j);
      end      
end
for i=1:numberNodesY
      for j=1:numberNodesX
        nodeCoordinates((i-1)*numberNodesX+j,2)=nodeCoordinatesY(i,1);
      end      
end
        
xx=nodeCoordinates(:,1);Lx=max(nodeCoordinatesX);
yy=nodeCoordinates(:,2);Ly=max(nodeCoordinatesY);
numberNodesX=size(nodeCoordinatesX,1);
numberNodesY=size(nodeCoordinatesY,1);
numberNodes=numberNodesX*numberNodesY;

for i=1:numberElementsX;
    for j=1:numberElementsY;
        ele_indx=i+numberElementsX*(j-1);
        elementNodes(ele_indx,1)=numberNodesX*2*(j-1)+i*2-1;
        elementNodes(ele_indx,5)=elementNodes(ele_indx,1)+1;
        elementNodes(ele_indx,2)=elementNodes(ele_indx,5)+1;
        elementNodes(ele_indx,8)=elementNodes(ele_indx,1)+numberNodesX;
        elementNodes(ele_indx,9)=elementNodes(ele_indx,8)+1;
        elementNodes(ele_indx,6)=elementNodes(ele_indx,9)+1;
        elementNodes(ele_indx,4)=elementNodes(ele_indx,8)+numberNodesX;
        elementNodes(ele_indx,7)=elementNodes(ele_indx,4)+1;
        elementNodes(ele_indx,3)=elementNodes(ele_indx,7)+1;
    end
end

% GDof: global number of degrees of freedom
GDof=5*numberNodes; 

% computation of the system stiffness matrix and force vector
%................................................................

% computation of stiffness matrix
% for Mindlin plate element

% K : stiffness matrix

K=zeros(GDof);

% Gauss quadrature for bending part
[gaussWeights,gaussLocations]=gaussQuadratureAdvanced('advanced');
 
% cycle for element
% cycle for element
for e=1:numberElements       
  % indice : nodal condofectivities for each element
  % elementDof: element degrees of freedom
  indice=elementNodes(e,:);           
  elementDof=[indice indice+numberNodes indice+2*numberNodes indice+3*numberNodes indice+4*numberNodes];    
  ndof=length(indice);
  L1=zeros(2*ndof,ndof*5);L1(1:ndof,3*ndof+1:4*ndof)=eye(9);L1(ndof+1:2*ndof,4*ndof+1:5*ndof)=eye(9);
  L2=zeros(2*ndof,ndof*5);L2(1:ndof,ndof+1:ndof*2)=eye(9);L2(ndof+1:2*ndof,2*ndof+1:ndof*3)=eye(9);
  F=zeros(ndof,ndof*5);F(1:ndof,1:ndof)=eye(9);
  
  % cycle for Gauss point
  for q=1:size(gaussWeights,1)                      
    GaussPoint=gaussLocations(q,:);                                                     
    xi=GaussPoint(1);
    eta=GaussPoint(2);

% shape functions and derivatives
    [shapeFunction,naturalDerivatives]=shapeFunctionQ9(xi,eta);

% Jacobian matrix, inverse of Jacobian 
% derivatives w.r.t. x,y    
    [Jacob,invJacobian,XYderivatives]=...
        Jacobian(nodeCoordinates(indice,:),naturalDerivatives);
    
% [S] matrix bending
    S=zeros(3,2*ndof);
    S(1,1:ndof)  = XYderivatives(:,1)';  
    S(2,ndof+1:2*ndof)= XYderivatives(:,2)';
    S(3,1:ndof)  = XYderivatives(:,2)';  
    S(3,ndof+1:2*ndof)= XYderivatives(:,1)';
 
% bending stiffness matrix 
    K(elementDof,elementDof)=K(elementDof,elementDof)+ ...
        L1'*S'*A*S*L1*gaussWeights(q)*det(Jacob)+...
        L1'*S'*B*S*L2*gaussWeights(q)*det(Jacob)+...
        L2'*S'*B*S*L1*gaussWeights(q)*det(Jacob)+...
        L2'*S'*D*S*L2*gaussWeights(q)*det(Jacob);
    
% [W] matrix bending (9*27)
    W=zeros(ndof,3*ndof);
    W(1:ndof,1:ndof)  = XYderivatives(:,1)*XYderivatives(:,1)';  
    W(1:ndof,1+ndof:2*ndof)  = XYderivatives(:,2)*XYderivatives(:,2)';
    W(1:ndof,1+2*ndof:3*ndof)  = XYderivatives(:,1)*XYderivatives(:,2)'...
                                +XYderivatives(:,2)*XYderivatives(:,1)';
                          
% [W1] matrix bending (27*9)                            
    W1=zeros(3*ndof,ndof);
    W1(1:ndof,1:ndof)=N_T(1,1)*eye(9);
    W1(1+ndof:2*ndof,1:ndof)=N_T(2,1)*eye(9);
    W1(1+2*ndof:3*ndof,1:ndof)=N_T(3,1)*eye(9);

% Thermal stiffness matrix
    K(elementDof,elementDof)=K(elementDof,elementDof)- ...
        F'*W*W1*F*gaussWeights(q)*det(Jacob);
    
% Gw matrix for VK stiffness
    
        w=zeros(ndof,1);
    
    Gw1=[XYderivatives(:,1)*XYderivatives(:,1)'*w];
    Gw2=[XYderivatives(:,2)*XYderivatives(:,2)'*w];
    Gw3=[XYderivatives(:,2)*XYderivatives(:,1)'*w]+[XYderivatives(:,1)*XYderivatives(:,2)'*w];
    Gw=[Gw1,Gw2,Gw3];
       
    end  % Gauss point
end   % element


% shear stiffness matrix

% Gauss quadrature for shear part
[gaussWeights,gaussLocations]=gaussQuadratureAdvanced('complete');
 
% cycle for element
% cycle for element
for e=1:numberElements       
  % indice : nodal condofectivities for each element
  % elementDof: element degrees of freedom
  indice=elementNodes(e,:);           
  elementDof=[indice indice+numberNodes indice+2*numberNodes indice+3*numberNodes indice+4*numberNodes];   
  ndof=length(indice);
  X=zeros(ndof,ndof*5);X(1:ndof,1:ndof)=eye(9);
  Y=L2;
  % cycle for Gauss point
  for q=1:size(gaussWeights,1)                      
    GaussPoint=gaussLocations(q,:);                                                     
    xi=GaussPoint(1);
    eta=GaussPoint(2);

% shape functions and derivatives
    [shapeFunction,naturalDerivatives]=shapeFunctionQ9(xi,eta);

% Jacobian matrix, inverse of Jacobian, 
% derivatives w.r.t. x,y    
    [Jacob,invJacobian,XYderivatives]=...
        Jacobian(nodeCoordinates(indice,:),naturalDerivatives);    
% [n] matrix shear
    n=zeros(2,ndof);
    n(1,1:ndof) = XYderivatives(:,1)';  
    n(2,1:ndof) = XYderivatives(:,2)';
% [H] matrix shear
    H=zeros(2,2*ndof);
    H(1,1:ndof) = shapeFunction'; 
    H(2,ndof+1:2*ndof) = shapeFunction';
    
    C=n*X+H*Y;

% stiffness matrix shear
    K(elementDof,elementDof)=K(elementDof,elementDof)+...
        C'*Ss*C*gaussWeights(q)*det(Jacob);  
  end  % gauss point
end    % element

stiffness=K;

%................................................................
% computation of  mass matrix 
% Ref.drv.#1

% mass : mass matrix 
mass=zeros(GDof);

% Gauss quadrature for bending part
[gaussWeights,gaussLocations]=gaussQuadratureAdvanced('advanced');
 
% cycle for element
for e=1:numberElements       
  % indice : nodal condofectivities for each element
  indice=elementNodes(e,:);
  elementDof=[indice indice+numberNodes indice+2*numberNodes indice+3*numberNodes indice+4*numberNodes];     
  ndof=length(indice);
  X=zeros(ndof,ndof*5);X(1:ndof,1:ndof)=eye(9);
  X1=zeros(ndof,ndof*5);X1(1:ndof,1+ndof:2*ndof)=eye(9);
  X2=zeros(ndof,ndof*5);X2(1:ndof,1+ndof*2:3*ndof)=eye(9);
  X3=zeros(ndof,ndof*5);X3(1:ndof,1+ndof*3:4*ndof)=eye(9);
  X4=zeros(ndof,ndof*5);X4(1:ndof,1+ndof*4:5*ndof)=eye(9);
  
  % cycle for Gauss point
  for q=1:size(gaussWeights,1)                      
    GaussPoint=gaussLocations(q,:);                                                       
    xi=GaussPoint(1);
    eta=GaussPoint(2);
    
% shape functions and derivatives
    [shapeFunction,naturalDerivatives]=shapeFunctionQ9(xi,eta);

% Jacobian matrix, inverse of Jacobian
% derivatives w.r.t. x,y    
    [Jacob,invJacobian,XYderivatives]=...
        Jacobian(nodeCoordinates(indice,:),naturalDerivatives);


% mass matrix
    N=shapeFunction';
    
    mass(elementDof,elementDof)=mass(elementDof,elementDof)+...
        (I0*(X'*N'*N*X + X3'*N'*N*X3 + X4'*N'*N*X4)+...
         I1*(X3'*N'*N*X1 + X4'*N'*N*X2)+...
         I2*(X1'*N'*N*X1 + X2'*N'*N*X2))*gaussWeights(q)*det(Jacob);
  end  % Gauss point
end    % element

%...................................................

% % % boundary conditions (ssss)
% fixedNodeW =find(xx==max(nodeCoordinates(:,1))|...
%                  xx==min(nodeCoordinates(:,1))|...
%                  yy==min(nodeCoordinates(:,2))|...
%                  yy==max(nodeCoordinates(:,2)));
% 
% fixedNodeTX =find(yy==max(nodeCoordinates(:,2))|...
%     yy==min(nodeCoordinates(:,2)));
% fixedNodeTY =find(xx==max(nodeCoordinates(:,1))|...
%     xx==min(nodeCoordinates(:,1)));
% fixedNodeU =find(xx==min(nodeCoordinates(:,1)));
% fixedNodeV =find(yy==min(nodeCoordinates(:,2)));
% prescribedDof=[fixedNodeW;fixedNodeTX+numberNodes;...
%       fixedNodeTY+2*numberNodes;...
%       fixedNodeU+3*numberNodes;fixedNodeV+4*numberNodes];
% activeDof=setdiff([1:GDof]',[prescribedDof]);

% % boundary conditions (cccc)
fixedNodeW =find(xx==max(nodeCoordinates(:,1))|...
                 xx==min(nodeCoordinates(:,1))|...
                 yy==min(nodeCoordinates(:,2))|...
                 yy==max(nodeCoordinates(:,2)));
fixedNodeTX =fixedNodeW;
fixedNodeTY =fixedNodeTX;
fixedNodeU =fixedNodeTX;
fixedNodeV =fixedNodeTX;
prescribedDof=[fixedNodeW;fixedNodeTX+numberNodes;...
      fixedNodeTY+2*numberNodes;...
      fixedNodeU+3*numberNodes;fixedNodeV+4*numberNodes];
activeDof=setdiff([1:GDof]',[prescribedDof]);

% solution
G=E_c/2.6;
% V : mode shape
% D : frequency
% 
E_m_v=2.0775e+011;
Rho_m_v=8166;
numberOfModes=2;
[V,D1] = eig(stiffness(activeDof,activeDof),...
    mass(activeDof,activeDof)); 
%D1 = diag(sqrt(D1)*pi*a^2/h*sqrt(Rho_m_v/E_m_v));
%shen verify
D1=diag(sqrt(D1)*a^2/pi^2*sqrt(h*Rho_m_v/(E_m_v*h^3/12/(1-poisson^2))));

[D1,ii] = sort(D1); ii = ii(1:numberOfModes); 
DD1(1:1,tk)=D1(1:1);
hoka=isreal(DD1);
    if hoka==0
       tkreal=tk;
       dT=dT(1:tk)';
       break
    end
   
end

DD1
DDD=DD1(1,:)'

