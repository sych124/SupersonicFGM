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
format short
tic

kt=[1];
k_interv = 0.025;
% Property of FGplate

Lratioi = [0.2/0.025];
Lratioe = [0.2/0.025];
Lratio_interv = 2.5;

a=0.2; Lx = a;
b=0.2; Ly = b;
thickness= a/Lratioi(1);h=thickness;
nt= [1000]; n_interv = 1;



time_interv = 1/640000;
time_end = 0.2;
time_range = 0:time_interv:time_end;

F_term = 0.0001;
F_amp = 1e6;

dT=[0]; dT_interv = 0.1;
Temperi = [400];
Tempere = [400];  %% Temperature of Tc
Temper_interv = 10;
temper_en = 1;
Temper_factor = 0;
delta_T = dT(1);

phi_air=0;               % angle of air flow
phi_air=phi_air*(pi/180);
C_phi=cos(phi_air);
S_phi=sin(phi_air);
lambda=600;          %Nondimension lambda = Lambda*Lx^4 / ( E_m*h^3 / ( 12*(1-nu^2) ) )
c_d=0.1;                                   % Mass_ratio coefficient (mu/Ma)     

MTS_en = 1;
disp('Predict Consumed Time =')
disp(1.9*size(dT(1):n_interv:dT(end)))
%%
for ROMtoMTS = 1:2
    MTS_en = ROMtoMTS-1;
    k_index = kt(1);
for Tm0 = dT(1):dT_interv:dT(end)
    step = nt(1);
    disp(k_index)
    disp('=')
for Temper = Temperi(1):Temper_interv:Tempere(1)

    
modeNumber=4;

T_ref=300;


Temper = T_ref + Tm0;
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
%  E_m_vec(2) = 105.7e9;
%  E_c_vec(2) = 320.24e9;
% nu_c_vec(2) = 0.26;
% nu_m_vec(2) = 0.2981;
% alpha_m_vec(2)=23e-6;
% alpha_c_vec(2)=7.4e-6;
% k_c = 10.4;
% k_m = 204;
% Rho_m_vec(2)= 4429;
% Rho_c_vec(2) = 3750;
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
%% Material Properties Pre-Setting Initializing
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
B = 0*BnMat_MTS; %Neutral surface Effect
D = h^3*DnMat_MTS;  

Ss = h*ShrCor_MTS*Ss_MTS*eye(2);

Nx = h*NnMMat_MTS;
Mx = MtE_MTS;

kappa=ShrCor_MTS; %ShrCor(n, E_mori, nu_mori, h0/h, z/h);
else
    A = h*AnMat_ROM;
    B = h^2*BnMat_ROM;
    B = 0*BnMat_ROM; %Neutral surface Effect
    D = h^3*DnMat_ROM;
    Ss = h*ShrCor_ROM*Ss_ROM*eye(2);
    
    Nx = h*NnMMat_ROM;
    Mx = MtE_ROM;
    
    kappa=ShrCor_ROM;
end

As = Ss;
N_T = Nx; N_delta_T = Nx;
M_T = Mx; M_delta_T = Mx;

I0 = MassI(1);
I1 = MassI(2);
I2 = MassI(3);


%% Element Connectivity
tk=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Gnerate Meshes %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    x_ele=5;                       % the number of elements for the X direction
    y_ele=5;                       % the number of elements for the Y direction
    x_node=x_ele*2+1;              % the number of nodes for the X direction
    y_node=y_ele*2+1;              % the number of nodes for the Y direction

    tot_ele=x_ele*y_ele;            % total number of elements
    tot_node=x_node*y_node;         % total number of nodes

%%%%%%%%%%%%%%%%%%%%%%%%%%%% define element type %%%%%%%%%%%%%%%%%%%%%%%
    ele_node=9;                     % the number of node per element
    node_dof=5;                     % D.O.F per node
    ele_dof=ele_node*node_dof;      % D.O.F per element
    tot_dof=tot_node*node_dof;      % total no of D.O.F 

    lx_e=Lx/x_ele;                  % x length of an element
    ly_e=Ly/y_ele;                  % y length of an element



%%%%%%%%%%%%%%%%%%%%%%  coordinate of the nodes  %%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:x_node
        for j=1:y_node
            n=i+(j-1)*x_node;
            coord(n,1)=(i-1)*(lx_e/2);
            coord(n,2)=(j-1)*(ly_e/2);
        end
    end

%%%%%%%%%%%%%%%%%%%%%% Boundary condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    BC_Simply_Supported
    %BC_Clamped

    bc_val=zeros(1,length(bc_dof));                 % Boundary values

%%%%%%%%%%%%%%%%%%% element connectivity %%%%%%%%%%%%%%%%%%%%%%%%%

    Connectivity                % node numbers for each element
% ex) node number of n'th node for m'th element is 28
%      --->  node_index(m,n)=28




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   PROCESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% Initialization of Matrix & Vector%%%%%%%%%%%%%%%%%%%%%%


    d=zeros(tot_dof,1);                     %displacement vector
    
%%%%%%%%%%%%%%%%%%%%%%% Gauss Quadrature %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Gauss_Quadrature
% 9pt integration for bending
% 4pt integration for shear


%%%%%%%%%%%%%%%%%%%%%%% Shape function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Shape_function;

% H(9,ordb)                 ordb=9
% dH_dxi(9,ordb)
% dH_deta (9,ordb)

% Hs(9,ords)                ords=4
% dHs_dxi(9,ords)            
% dHs_deta (9,ords)


%%%%%%%%%%%%%%%%%%%%% Linear Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Linear_Matrix;

    
    Thermal_Load;
% You can get  K_elastic & K_thermal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Apply Boundary Condition %%%%%%%%%%%%%%%%%%%%

    bc_length=length(bc_dof);

    for i=1:bc_length
        for j=1:tot_dof
            K_elastic(bc_dof(i),j)=0;
            K_thermal(bc_dof(i),j)=0;
            K_elastic(j,bc_dof(i))=0;
            K_thermal(j,bc_dof(i))=0;
        end
    
        K_elastic(bc_dof(i),bc_dof(i))=1;
        K_thermal(bc_dof(i),bc_dof(i))=1;
    
    end
%%     Calculation of Critical Buckling Temperature (Linear Analsis)   

    [V,d] = eig(K_elastic,K_thermal) ;
 
    [cr_T, k] = sort(diag(d)) ; 

    V=V(:,k) ;

 
    T= real(cr_T)  ;
  
    Temp_cr = T(min(find(T>1)))
 
     
%% FEM  Calculations for  1st & 2nd Nonlinear Stiffness matrix  
  
    %-------------------------initial displacement for N-R method ------------                           

    EV=real(V(:,min(find(T>1))));
 
    w_vec = EV(tot_node*4+1:tot_node*5,1) ;
    
    z= zeros(y_node,x_node);
  
 %-------Mode Shape for thermal post buckling------------------------------    
    for r=1:x_node 
        bb = w_vec(1+(r-1)*y_node : y_node+(r-1)*y_node)/h ;
        z(r,:) = bb' ;
    end

    x= [0:(Lx/(x_node-1)):Lx ];
    y= [0:(Ly/(y_node-1)):Ly ];

    max_w = max(abs(w_vec)) ;
    D_m=Em_eff*h^3/(12*(1-(mean(poissn_ROM))^2));          %flexural rigidity of metal
    Lambda=lambda/(Lx^3/D_m);     % converge nondimensionalized dynamic pressere to the actual value of dynamic pressure
    
    scale = -1.0;                    %scale factor for initial condition
    
    disp_s = scale*h*abs(EV)/max_w ;    % the initial estimated deflection 
    
    
    itr=1;
    error_norm=1;
   
    while(error_norm>=1.0e-7)
        
        itr=itr
                
        Static_Nonlinear_Matrix;   %run Static_Nonlinear_Matrix.m
                                   %and get N1_s & N2_s
        Tangential_stiffness=zeros(tot_dof,tot_dof);
        
        Tangential_stiffness=K_elastic - delta_T*K_thermal + Lambda*A_F + N1_s + N2_s;
      
        Load_imbalance=delta_T*F_thermal-(K_elastic - delta_T*K_thermal + Lambda*A_F + 1/2*N1_s + 1/3*N2_s)*disp_s;

        %------------------------- Apply Boundary Condition --------------

        bc_length=length(bc_dof);

        for j=1:bc_length
            for k=1:tot_dof
                Tangential_stiffness(bc_dof(j),k)=0;
            end
            Tangential_stiffness(bc_dof(j),bc_dof(j))=1;
            Load_imbalance(bc_dof(j),1)=bc_val(j);
        end
        
      
        delta_disp=inv(Tangential_stiffness)*Load_imbalance;
        

        %Convergence check
        
        if(itr==1)
            error_norm=1.0
        else
            error_norm=norm(delta_disp(4*tot_node+1:5*tot_node))/h
        end
                
        disp_s=disp_s+delta_disp;
        itr=itr+1;
  
        if(itr==40)
            break            
        end
    end

    max_disp_s = max(abs(disp_s(4*tot_node+1:5*tot_node)));
    max_ratio_s=max_disp_s/h;
 
% -------------------------plot displacement shape----------------------------
    z= zeros(y_node,x_node);
  
    w_result=disp_s(4*tot_node+1:5*tot_node);
    
   for r=1:x_node 
        for q=1:y_node
            z(q,r)=w_result((q-1)*x_node+r);
        end
   end
    
    x= [0:(Lx/(x_node-1)):Lx ];
    y= [0:(Ly/(y_node-1)):Ly ];
 

%% solution of Time-Domain
NDof = size(activeDof); NDof = NDof(1,1);

init_def = zeros(1,NDof);
init_vel = zeros(1,NDof);
F_init = [ F_amp*zeros(44,1);F_amp*ones(11,1);F_amp*zeros(433,1);F_amp*ones(11,1);F_amp*zeros(0,1)]; % foc SS  [ zeros(NDof/3,1);F_amp*ones(NDof/3,1);zeros(NDof/3,1)];  for CC %
F = [F_init.*ones(NDof,int64(F_term/time_interv)) zeros(NDof,int64(time_end/time_interv)-int64(F_term/time_interv))];
deform_time = newmark(mass(activeDof,activeDof), 0*mass(activeDof,activeDof), (stiffness(activeDof,activeDof)-K_T(activeDof,activeDof)), F,time_range,init_def, init_vel);

dml_deform_time = deform_time/h;

figure(1);
plot(time_range(1:2500),dml_deform_time(1:2500,int64(end-121)),'LineWidth',2);
grid on;
hold on;
%%
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


[V,D2] = eig(stiffness(activeDof,activeDof),-K_T(activeDof,activeDof));
D2 = diag(D2);
D2 = sort(D2);

% Vibration
Da = sort(Da);
disp(Da(1))


end

Da_k(:,int64((Tm0-dT(1))/dT_interv)+1) = sort(Da);
freq_k(:,int64((Tm0-dT(1))/dT_interv)+1) = sort(freq);

end


sample_interv = time_interv;

time = time_range(1001:end-1);
Pren = transpose(deform_time(1001:end,end-121));

Fs = 1/sample_interv;
[L,s] = size(Pren');

Y = fft(Pren);
f = Fs*(0:(L/2))/L;
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
%% Find Max
maxP1 = max(P1);
minP1 = min(P1);
indexOfFirstMax = find(P1 == maxP1, 1, 'first');
maxY = P1(indexOfFirstMax);
maxX = f(indexOfFirstMax);
maxXYFs = [maxX, maxY,Fs];
%% Ploting on decade scale
%{
figure;
plot(f,P1,'LineWidth',2); %%for decade scale plot
grid on;
axis([maxX/5 500 minP1 maxY*1.15]);
tstr = [fp ,' FFT Signal',' (',num2str(TimeS),'sec ~ ',num2str(TimeE),'sec)'];
title(tstr,'FontSize',25);
xlabel('f (Hz)');
ylabel('Y(t)');
str = [' <- max (' num2str(maxX) ',' num2str(maxY) ')']; 
text(maxX,maxY,str,'FontSize',18);
%}
%% Ploting on log scale
figure(2);
loglog(f,P1,'-','LineWidth',2);
axis([10e1 20e4 10e-7 10]);
grid on;
xlabel('f (Hz)');
ylabel('Y(t)');
str = [' <- max (' num2str(maxX) ',' num2str(maxY) ')']; 
text(maxX,maxY,str,'FontSize',18);
hold on;

DD1_RM(:,ROMtoMTS) = D1;
DD2_RM(:,ROMtoMTS) = D2;
freq_RM(:,ROMtoMTS) = sort(freq);
DDa_RM(:,ROMtoMTS) = sort(Da);

end

%plot(dT(1):dT_interv:dT(end),freq_k(1,:));

consumed_time = toc