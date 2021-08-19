%................................................................
 
% 1st Paper Project in 2015 
% Vibaratino analysis of functionally graded plates in thermal effects
% Written by Kang

% Concept

% 1. Nonlinear/uniform temperature distribution
% 2. Voigt-rule 
% 3. S-FGM
% 4. Neutral plane
% 5. von-Karman nonlinear term
% 6. Based on first-order shear deformation plate theory


% clear memory
clear all;colordef white;
clc
format short
tic
result = zeros(1,length(0:200));

kt=[0.55 1];
k_interv = 0.05;
% Property of FGplate

Lratioi = [40];
Lratioe = [40];
Lratio_interv = 2.5;

a=0.4; % a=Lx
b=0.4; % b=Ly
thickness= a/Lratioi(1);h=thickness;
nt= [10000]; n_interv = 1;
step = nt(1);


dT=[0];
Temperi = [350];
Tempere = [350];  %% Temperature of Tc
Tm0 = 50; %% Temperature difference T_ref and Tm
Temper_interv = 0.01;
temper_en = 1;
Temper_factor = 0;
ShrTemp = 0;
dT=[Temperi(1)];

MTS_en = 1;
disp('Predict Consumed Time =')
disp(19*size(nt(1):n_interv:nt(end)))
for ROMtoMTS = 1:1
    MTS_en = ROMtoMTS-1;
    
for k_index = kt(1):k_interv:kt(end)
    disp(k_index)
    disp('=')
for Temper = Temperi(1):Temper_interv:Tempere(1)

modeNumber=4;

T_ref=300;
T_delta = Temper - T_ref - Temper_factor;
T0 = T_ref;
Tm = T_ref+Tm0;
Tc = Temper;

prach = h;

I=h^3/12;

[E_c_vec, alpha_c_vec,Rho_c_vec,  k_c_vec, nu_c_vec] = MatTemp('SiN');
[E_m_vec, alpha_m_vec, Rho_m_vec, k_m_vec, nu_m_vec] = MatTemp('SUS');
if temper_en == 0
%  Em_eff= 67e9;
%  Ec_eff= 302e9*0.4+0.6*67e9;
%  nu_m = 0.33;
% nu_c = 0.17*0.4+0.33*0.6;
% alpha_m=23e-6;
% alpha_c=7.4e-6;
% k_c = 10.4;
% k_m = 204;
%  rho_m = 2700;
%  rho_c = 3200*0.4+2700*0.6;
 E_m_vec(2) = 70e9;
 E_c_vec(2) = 380e9;
nu_c_vec(2) = 0.3;
nu_m_vec(2) = 0.3;
alpha_m_vec(2)=23e-6;
alpha_c_vec(2)=7.4e-6;
k_c = 10.4;
k_m = 204;
Rho_m_vec(2)= 2702;
Rho_c_vec(2) = 3800;
end
%% Material Properties

Area=1;
%
P = 0; % uniform pressure
% constitutive matrix 
% E=[1 poisson 0;poisson 1 0;0 0 (1-poisson)/2]/(1-poisson^2);
 %% Pre-Setting Initializing
interv = 1/(step);
z = -h/2 : h/(step) : prach/2;% Practical Z direction
dmlZ = z/h; % dimensionless Z (Z_d),
midZ = dmlZ;%dmlZ(1:end-1) + interv/2; %mid-line of each step
dmlV = dmlZ;%[-0.5, midZ(2:end-1), 0.5]; %dimensionless Fraction
Vc = (dmlV+1/2).^k_index; Vm = 1-Vc; % Volume Fraction

%% Initial contitions
%% Temperature Distribution
for temp = 1:step+1
[Kc_eff, Km_eff, K_ROM] = MatProp(k_c_vec, k_m_vec,T_ref+T_delta,k_index,dmlV,temper_en);

Kcm = Kc_eff-Km_eff;
P = z(temp)/h+0.5;
C=1-Kcm/((k_index+1)*Km_eff)+Kcm^2/((2*k_index+1)*Km_eff^2)-Kcm^3/((3*k_index+1)*Km_eff^3)+Kcm^4/((4*k_index+1)*Km_eff^4)...
        -Kcm^5/((5*k_index+1)*Km_eff^5);
Rz=(P-Kcm/((k_index+1)*Km_eff)*P^(k_index+1)+Kcm^2/((2*k_index+1)*Km_eff^2)*P^(2*k_index+1)...
        -Kcm^3/((3*k_index+1)*Km_eff^3)*P^(3*k_index+1)+Kcm^4/((4*k_index+1)*Km_eff^4)*P^(4*k_index+1)...
        -Kcm^5/((5*k_index+1)*Km_eff^5)*P^(5*k_index+1))/C;
    
Tz(temp) = Tm+(Tc-Tm)*Rz;

end

%% Temperature-dependent Material Properties
for temp = 1 : step+1
[Ec_eff, Em_eff, Ent_ROM] = MatProp(E_c_vec,E_m_vec,Tz(temp),k_index,dmlV,temper_en);
[poissc_eff, poissm_eff, poissnt_ROM] = MatProp(nu_c_vec, nu_m_vec,Tz(temp),k_index,dmlV,temper_en);
[rhoc_eff, rhom_eff, rhont_ROM] = MatProp(Rho_c_vec, Rho_m_vec,Tz(temp),k_index,dmlV,temper_en);
[alphc_eff, alphm_eff, alpht_ROM] = MatProp(alpha_c_vec, alpha_m_vec,Tz(temp),k_index,dmlV,temper_en);

En_ROM(temp) = Ent_ROM(temp);
poissn_ROM(temp) = poissnt_ROM(temp);
rhon_ROM(temp) = rhont_ROM(temp);
alph_ROM(temp) = alpht_ROM(temp);

%poissn_ROM = 0*poissn_ROM + 0.28; %% Shen Wang

end
%% Computated Parameter 1 - Voigt Model
Epn_ROM = En_ROM;
Gpn_ROM = En_ROM./(2*(1.+poissn_ROM));
Gc_eff = Gpn_ROM(end); Gm_eff = Gpn_ROM(1);

z0n_ROM =  (midZ*transpose(Epn_ROM))/(sum(Epn_ROM));
z0_ROM = z0n_ROM .* h;
ShrCor_ROM =  ShrCor(step, En_ROM, poissn_ROM, z0n_ROM,dmlZ);
if ShrTemp == 0
    [Ect_eff, Emt_eff, Ent_ROM] = MatProp(E_c_vec,E_m_vec, T_ref,k_index,dmlV,temper_en);
    [poissct_eff, poissmt_eff, poissnt_ROM] = MatProp(nu_c_vec, nu_m_vec,T_ref,k_index,dmlV,temper_en);
    ShrCor_ROM = ShrCor(step, Ent_ROM, poissnt_ROM, z0n_ROM,dmlZ);
end

%% Cmoputed parameter 2 - Mori-Tanaka Model
Bc_eff = Ec_eff/(3*(1-2*poissc_eff));
Bm_eff = Em_eff/(3*(1-2*poissm_eff));

f1 = Gm_eff*(9*Bm_eff + 8*Gm_eff)/(6*(Bm_eff + 2*Gm_eff));

Bpn_MTS = (Vc./(1+Vm.*((3*Bc_eff-3*Bm_eff)/(3*Bm_eff+4*Gm_eff)))) ...
    .*(Bc_eff-Bm_eff) + Bm_eff;
Gpn_MTS = (Vc./(1+Vm.*((Gc_eff-Gm_eff)/(Gm_eff+f1)))) ...
    .*(Gc_eff-Gm_eff) + Gm_eff;

alph_MTS = (1./Bpn_MTS - 1/Bm_eff)/(1/Bc_eff - 1/Bm_eff) ...
        .*(alphc_eff-alphm_eff) + alphm_eff;

for temp = 1 : step+1
    Epn_MTS(temp) = 9*Bpn_MTS(temp)*Gpn_MTS(temp)/(3*Bpn_MTS(temp)+Gpn_MTS(temp));
    poissn_MTS(temp) = (3*Bpn_MTS(temp)-2*Gpn_MTS(temp))/(2*(3*Bpn_MTS(temp)+Gpn_MTS(temp)));
end
z0n_MTS = (midZ*transpose(Epn_MTS))/(sum(Epn_MTS));
z0_MTS = z0n_MTS .* h;

%poissn_MTS = 0*poissn_MTS + 0.28; %% Shen Wang
ShrCor_MTS = ShrCor(step, Epn_MTS,poissn_MTS, z0n_MTS,dmlZ);
%% Computated parameter 2 - nondimensional ABD NM with Voigt
Emat=[1 0.3 0;0.3 1 0;0 0 (1-0.3)/2]; %[1 nu 0; nu 1 0; 0 0 (1-nu)/2]  %  (Does not used)

AnE_ROM = [sum(interv.*Epn_ROM); sum(interv.*Gpn_ROM)]; %  (Does not used)
AnE1_ROM = sum(interv.*Epn_ROM./(1-poissn_ROM.^2)); 
AnE2_ROM = interv.*(poissn_ROM./(1-poissn_ROM.^2))*transpose(Epn_ROM); 
AnE3_ROM = interv*sum(Gpn_ROM);
AnMat_ROM = [AnE1_ROM, AnE2_ROM, 0; AnE2_ROM, AnE1_ROM, 0; 0, 0, AnE3_ROM];

BnE_ROM = interv.*(midZ - z0n_ROM)*transpose([Epn_ROM;Gpn_ROM]);  %  (Does not used)
BnE1_ROM = interv.*(midZ - z0n_ROM)*transpose(Epn_ROM./(1-poissn_ROM.^2));
BnE2_ROM = interv.*(midZ - z0n_ROM).*(poissn_ROM./(1-poissn_ROM.^2))*transpose(Epn_ROM); 
BnE3_ROM = interv.*(midZ - z0n_ROM)*transpose(Gpn_ROM);
BnMat_ROM = [BnE1_ROM, BnE2_ROM, 0; BnE2_ROM, BnE1_ROM, 0; 0, 0, BnE3_ROM];

DnE_ROM = interv.*((midZ     - z0n_ROM).^2)*transpose([Epn_ROM;Gpn_ROM]);  %  (Does not used)
DnE1_ROM = interv.*((midZ     - z0n_ROM).^2)* transpose(Epn_ROM./(1-poissn_ROM.^2)); 
DnE2_ROM = interv.*((midZ     - z0n_ROM).^2).*(poissn_ROM./(1-poissn_ROM.^2))*transpose(Epn_ROM); 
DnE3_ROM = interv.*((midZ     - z0n_ROM).^2) * transpose(Gpn_ROM);
DnMat_ROM = [DnE1_ROM, DnE2_ROM, 0; DnE2_ROM, DnE1_ROM, 0; 0, 0, DnE3_ROM];

Ss_ROM = sum(interv.*Epn_ROM./(2*(1+poissn_ROM))); 

NtnE_ROM = interv*(Tz-T0).*alph_ROM*transpose(Epn_ROM); 
Ntn1E_ROM = interv*(Tz-T0).*alph_ROM*transpose((poissn_ROM./(1-poissn_ROM.^2)).*Epn_ROM); 
Ntn2E_ROM = interv*(Tz-T0).*alph_ROM*transpose((poissn_ROM./(1-poissn_ROM.^2)).*Epn_ROM); 
Ntn3E_ROM = interv*(Tz-T0).*alph_ROM*transpose(Gpn_ROM); 
NnMMat_ROM = [Ntn1E_ROM; Ntn2E_ROM; Ntn3E_ROM;];

MtnE_ROM = interv*(Tz-T0).*(midZ-z0n_ROM)*diag(alph_ROM)*transpose(Epn_ROM);  
%% Computated parameter 3 - dimensional ABD NM with Voigt
AE_ROM = h*AnE_ROM; 
BE_ROM = h^2*BnE_ROM;
DE_ROM = h^3*DnE_ROM;
GE_ROM = h*ShrCor_ROM*AnE_ROM(2);
NtE_ROM = h*NtnE_ROM;
MtE_ROM = h^2*MtnE_ROM;
SE_ROM = h*interv.*(sum(alph_ROM).*(Tm-T0)).^2;
%% Computed parameter 4 - transform to FEM  (Does not used)
% AMat = AE_ROM(1,1)*Emat; As(1,1) = AMat(3,3); As(2,2) = AMat(3,3);
% BMat = BE_ROM(1,1)*Emat;
% DMat = DE_ROM(1,1)*Emat;
% 
% N_delta_T=zeros(3,1);
% M_delta_T=zeros(3,1);
% 
% N_delta_T(1,1)=NtE_ROM;
% N_delta_T(2,1)=N_delta_T(1,1);
% 
% %M_delta_T(1,1)=h^2/(1-nu) * (1/2*(E_c-E_m)*(alpha_c-alpha_m)*(k_index/((k_index+1)*(2*k_index+1))) + (alpha_m*(E_c-E_m)+E_m*(alpha_c-alpha_m))*k_index/(2*(k_index+2)*(k_index+1)));
% M_delta_T(1,1)=MtE_ROM;
% M_delta_T(2,1)=M_delta_T(1,1);

MassIn = [sum(interv.*rhon_ROM);
    interv.*(midZ - z0n_ROM)*transpose(rhon_ROM);
    interv.*((midZ - z0n_ROM).^2)*transpose(rhon_ROM)];
MassI = [h*MassIn(1);h^2*MassIn(2);h^3*MassIn(3)]; 


%% Computated parameter 2 - nondimensional ABD NM with MTS
Emat=[1 0.3 0;0.3 1 0;0 0 (1-0.3)/2];

AnE_MTS = [sum(interv.*Epn_MTS); sum(interv.*Gpn_MTS)]; %  (Does not used)
AnE1_MTS = sum(interv.*Epn_MTS./(1-poissn_MTS.^2)); 
AnE2_MTS = interv.*(poissn_MTS./(1-poissn_MTS.^2))*transpose(Epn_MTS); 
AnE3_MTS = interv*sum(Gpn_MTS);
AnMat_MTS = [AnE1_MTS, AnE2_MTS, 0; AnE2_MTS, AnE1_MTS, 0; 0, 0, AnE3_MTS];

BnE_MTS = interv.*(midZ - z0n_MTS)*transpose([Epn_MTS;Gpn_MTS]);  %  (Does not used)
BnE1_MTS = interv.*(midZ - z0n_MTS)*transpose(Epn_MTS./(1-poissn_MTS.^2));
BnE2_MTS = interv.*(midZ - z0n_MTS).*(poissn_MTS./(1-poissn_MTS.^2))*transpose(Epn_MTS); 
BnE3_MTS = interv.*(midZ - z0n_MTS)*transpose(Gpn_MTS);
BnMat_MTS = [BnE1_MTS, BnE2_MTS, 0; BnE2_MTS, BnE1_MTS, 0; 0, 0, BnE3_MTS];

DnE_MTS = interv.*((midZ     - z0n_MTS).^2)*transpose([Epn_MTS;Gpn_MTS]);  %  (Does not used)
DnE1_MTS = interv.*((midZ     - z0n_MTS).^2)* transpose(Epn_MTS./(1-poissn_MTS.^2)); 
DnE2_MTS = interv.*((midZ     - z0n_MTS).^2).*(poissn_MTS./(1-poissn_MTS.^2))*transpose(Epn_MTS); 
DnE3_MTS = interv.*((midZ     - z0n_MTS).^2) * transpose(Gpn_MTS);
DnMat_MTS = [DnE1_MTS, DnE2_MTS, 0; DnE2_MTS, DnE1_MTS, 0; 0, 0, DnE3_MTS];

Ss_MTS = sum(interv.*Epn_MTS./(2*(1+poissn_MTS))); 

NtnE_MTS = interv*(Tz-T0).*alph_MTS*transpose(Epn_MTS);  %  (Does not used)
Ntn1E_MTS = interv*(Tz-T0).*alph_MTS*transpose((poissn_MTS./(1-poissn_MTS.^2)).*Epn_MTS); 
Ntn2E_MTS = interv*(Tz-T0).*alph_MTS*transpose((poissn_MTS./(1-poissn_MTS.^2)).*Epn_MTS); 
Ntn3E_MTS = interv*(Tz-T0).*alph_MTS*transpose(Gpn_MTS); 
NnMMat_MTS = [Ntn1E_MTS; Ntn2E_MTS; Ntn3E_MTS;];

MtnE_MTS = interv*(Tz-T0).*(midZ-z0n_MTS)*diag(alph_MTS)*transpose(Epn_MTS);  %  (Does not used)

%% Computated parameter 3 - dimensional ABD NM with MTS (does not used)
AE_MTS = h*AnE_MTS; 
BE_MTS = h^2*BnE_MTS;
DE_MTS = h^3*DnE_MTS; 
GE_MTS = h*ShrCor_MTS*AnE_MTS(2);
NtE_MTS = h*NtnE_MTS;
MtE_MTS = h^2*MtnE_MTS;
SE_MTS = h*interv.*(sum(alph_MTS).*(Tm-T0)).^2;
%% Computed parameter 4 - transform to FEM Mori-Tanaka
% AMat = AE_MTS(1,1)*Emat; As(1,1) = AMat(3,3); As(2,2) = AMat(3,3);
% BMat = BE_MTS(1,1)*Emat;
% DMat = DE_MTS(1,1)*Emat;

% N_delta_T=zeros(3,1);
% M_delta_T=zeros(3,1);
% 
% N_delta_T(1,1)=NtE_MTS;
% N_delta_T(2,1)=N_delta_T(1,1);
% 
% %M_delta_T(1,1)=h^2/(1-nu) * (1/2*(E_c-E_m)*(alpha_c-alpha_m)*(k_index/((k_index+1)*(2*k_index+1))) + (alpha_m*(E_c-E_m)+E_m*(alpha_c-alpha_m))*k_index/(2*(k_index+2)*(k_index+1)));
% M_delta_T(1,1)=MtE_MTS;
% M_delta_T(2,1)=M_delta_T(1,1);

MassIn = [sum(interv.*rhon_ROM);
    interv.*(midZ - z0n_ROM)*transpose(rhon_ROM);
    interv.*((midZ - z0n_ROM).^2)*transpose(rhon_ROM)];
MassI = [h*MassIn(1);h^2*MassIn(2);h^3*MassIn(3)]; 


%% Pre-Processing for FEM input
if MTS_en == 1
A = h*AnMat_MTS;
B = h^2*BnMat_MTS;
D = h^3*DnMat_MTS;

Ss = h*ShrCor_MTS*Ss_MTS*eye(2);

Nx = h*NnMMat_MTS;
Mx = MtE_MTS;

kapa=ShrCor_MTS; %ShrCor(n, E_mori, nu_mori, h0/h, z/h);
else
    A = h*AnMat_ROM;
    B = h^2*BnMat_ROM; 
    D = h^3*DnMat_ROM;
    Ss = h*ShrCor_ROM*Ss_ROM*eye(2);
    
    Nx = h*NnMMat_ROM;
    Mx = MtE_ROM;
    
    kapa=ShrCor_ROM;
end
N_T = Nx;

I0 = MassI(1);
I1 = MassI(2);
I2 = MassI(3);
%% Element Connectivity
for tk=1:size(dT,2);

%Mesh generation
% Generation of coordinates(9 node)
Lx=a;Ly=b;
numberElementsX=5;
numberElementsY=5;
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

K_T=zeros(GDof);
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
    
% [W1] matrix bending (9*27)
    W1=zeros(ndof,3*ndof);
    W1(1:ndof,1:ndof)  = XYderivatives(:,1)*XYderivatives(:,1)';  
    W1(1:ndof,1+ndof:2*ndof)  = XYderivatives(:,2)*XYderivatives(:,2)';
    W1(1:ndof,1+2*ndof:3*ndof)  = XYderivatives(:,1)*XYderivatives(:,2)'...
                                +XYderivatives(:,2)*XYderivatives(:,1)';
                          
% [W2] matrix bending (27*9)                            
    W2=zeros(3*ndof,ndof);
    W2(1:ndof,1:ndof)=N_T(1,1)*eye(9);
    W2(1+ndof:2*ndof,1:ndof)=N_T(2,1)*eye(9);
    W2(1+2*ndof:3*ndof,1:ndof)=N_T(3,1)*eye(9);
 
   
% Thermal stiffness matrix
    K_T(elementDof,elementDof)=K_T(elementDof,elementDof)- ...
        F'*W1*W2*F*gaussWeights(q)*det(Jacob);
    
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
fixedNodeW =find(xx==max(nodeCoordinates(:,1))|...
                 xx==min(nodeCoordinates(:,1))|...
                 yy==min(nodeCoordinates(:,2))|...
                 yy==max(nodeCoordinates(:,2)));

fixedNodeTX =find(yy==max(nodeCoordinates(:,2))|...
    yy==min(nodeCoordinates(:,2)));
fixedNodeTY =find(xx==max(nodeCoordinates(:,1))|...
    xx==min(nodeCoordinates(:,1)));
fixedNodeU =find(xx==min(nodeCoordinates(:,1)));
fixedNodeV =find(yy==min(nodeCoordinates(:,2)));
prescribedDof=[fixedNodeW;fixedNodeTX+numberNodes;...
      fixedNodeTY+2*numberNodes;...
      fixedNodeU+3*numberNodes;fixedNodeV+4*numberNodes];
activeDof=setdiff([1:GDof]',[prescribedDof]);

% % boundary conditions (cccc)
% fixedNodeW =find(xx==max(nodeCoordinates(:,1))|...
%                  xx==min(nodeCoordinates(:,1))|...
%                  yy==min(nodeCoordinates(:,2))|...
%                  yy==max(nodeCoordinates(:,2)));
% fixedNodeTX =fixedNodeW;
% fixedNodeTY =fixedNodeTX;
% fixedNodeU =fixedNodeTX;
% fixedNodeV =fixedNodeTX;
% prescribedDof=[fixedNodeW;fixedNodeTX+numberNodes;...
%       fixedNodeTY+2*numberNodes;...
%       fixedNodeU+3*numberNodes;fixedNodeV+4*numberNodes];
% activeDof=setdiff([1:GDof]',[prescribedDof]);

% solution

% V : mode shape
% D : frequency
% 


[rhoc_300, rhom_300, rhon300_ROM] = MatProp(Rho_c_vec, Rho_m_vec,300,k_index,dmlV,temper_en);
[Ec300f, Em_300, En300_ROM] = MatProp(E_c_vec,E_m_vec,300,k_index,dmlV,temper_en);

numberOfModes=1;
[V,D1] = eig(stiffness(activeDof,activeDof)+K_T(activeDof,activeDof),...
    mass(activeDof,activeDof)); 
freq = diag(sqrt(D1))/(2*pi);
Dn = diag(sqrt(D1))*a^2/pi^2*sqrt(I0/D(1));
if MTS_en ==0
    Da = diag(sqrt(D1))*(a^2/h)*sqrt((1-mean(poissn_ROM)^2)*rhom_eff/Em_eff); %Em_300 ??
else
    Da = diag(sqrt(D1))*(a^2/h)*sqrt((1-mean(poissn_MTS)^2)*rhom_eff/Em_eff); %Em_300 ??
end
D1 = diag(sqrt(D1)*h*sqrt(rhoc_eff/Ec_eff)); Dn = D1;
[D1,ii] = sort(D1); ii = ii(1:numberOfModes); 
DD1(1:1,tk)=D1(1:1);
hoka=isreal(DD1);
    if hoka==0
       tkreal=tk;
       dT=dT(1:tk)';
       break
    end
end

[V,D2] = eig(stiffness(activeDof,activeDof),-K_T(activeDof,activeDof));
D2 = diag(D2);
D2 = sort(D2);

% Vibration
Da = sort(Da);
disp(Da(1))


end

Da_k(:,int64((k_index-kt(1))/k_interv)+1) = sort(Da);
freq_k(:,int64((k_index-kt(1))/k_interv)+1) = sort(freq);

end
DD1_RM(:,ROMtoMTS) = D1;
DD2_RM(:,ROMtoMTS) = D2;
freq_RM(:,ROMtoMTS) = sort(freq);
DDa_RM(:,ROMtoMTS) = sort(Da);

end

consumed_time = toc