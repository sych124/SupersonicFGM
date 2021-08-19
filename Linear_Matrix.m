%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Linear_Matrix.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% You can get " elastic_stiffness & thermal_stiffness
%         & aerodynamic stiffness % aerodynamic damping matrix " from this


K_elastic=zeros(tot_dof,tot_dof);        % elastic stiffness matrix
K_thermal=zeros(tot_dof,tot_dof);        % thermal stiffness matrix 
A_F=zeros(tot_dof,tot_dof);              % aerodynamic stiffness matrix
A_D=zeros(tot_dof,tot_dof);              % aerodynamic damping matrix

for ne=1:tot_ele                                       

    for i=1:ele_node
        X_coord(i)=coord(node_indx(ne,i),1);
        Y_coord(i)=coord(node_indx(ne,i),2);
    end
%-------------------------------------------------------------------  
%                           3rd order integration
%-------------------------------------------------------------------
    kb=zeros(ele_dof,ele_dof);
    k_th=zeros(ele_dof,ele_dof);
    a_f=zeros(ele_dof,ele_dof);
    a_d=zeros(ele_dof,ele_dof);

    dx_dxi=zeros(1,ordb);
    dx_deta=zeros(1,ordb);
    dy_dxi=zeros(1,ordb);
    dy_deta=zeros(1,ordb);       %initialize summation terms
    
    detJ_b=zeros(1,ordb);            %initialize summation terms
    
    for i=1:ordb 
        for j=1:ele_node
           dx_dxi(1,i)= dx_dxi(1,i)+X_coord(j)*dH_dxi(j,i);
           dy_deta(1,i)= dy_deta(1,i)+Y_coord(j)*dH_deta(j,i);
           dx_deta(1,i)= dx_deta(1,i)+X_coord(j)*dH_deta(j,i);
           dy_dxi(1,i)= dy_dxi(1,i)+Y_coord(j)*dH_dxi(j,i);
        end
        
        detJ_b(1,i)=dx_dxi(1,i)*dy_deta(i)-dx_deta(i)*dy_dxi(i);    % Get Jacobian %
%    -------------------  banding stiffness matrix-------------------
        Bm=zeros(3,ele_dof);            
        Bb=zeros(3,ele_dof);            
        
        R=zeros(3,4);
        B_1=zeros(4,ele_dof);
        B_2=zeros(4,ele_dof);
        
        R(1,1)=dy_deta(1,i);
        R(1,2)=-dy_dxi(1,i);
        R(2,3)=-dx_deta(1,i);
        R(2,4)=dx_dxi(1,i);
        R(3,1)=R(2,3);
        R(3,2)=R(2,4);
        R(3,3)=R(1,1);
        R(3,4)=R(1,2);
        
        for j=1:ele_node
            B_1(1,2*ele_node+j)=dH_dxi(j,i);
            B_1(2,2*ele_node+j)=dH_deta(j,i);
            B_1(3,3*ele_node+j)=dH_dxi(j,i);
            B_1(4,3*ele_node+j)=dH_deta(j,i);
            
            B_2(1,j)=dH_dxi(j,i);
            B_2(2,j)=dH_deta(j,i);
            B_2(3,ele_node+j)=dH_dxi(j,i);
            B_2(4,ele_node+j)=dH_deta(j,i);
        end
        
        Bm=(R*B_1)/detJ_b(1,i);
        Bb=(R*B_2)/detJ_b(1,i);
        
        
        kb=kb+GaussWt_b(1,i)*detJ_b(1,i)*(Bm'*A*Bm + Bb'*B*Bm + Bm'*B*Bb + Bb'*D*Bb);
        
%    -------------------  Thermal stiffness matrix -------------------        
        Bth=zeros(2,ele_dof);            
        N_th=zeros(2,2);
      
        N_th=[ N_delta_T(1,1),N_delta_T(3,1);
               N_delta_T(3,1),N_delta_T(2,1)];
        
        for j=1:ele_node
            Bth(1,4*ele_node+j)=(dH_dxi(j,i)*dy_deta(1,i)-dH_deta(j,i)*dy_dxi(1,i))/detJ_b(1,i);
            Bth(2,4*ele_node+j)=(-dH_dxi(j,i)*dx_deta(1,i)+dH_deta(j,i)*dx_dxi(1,i))/detJ_b(1,i);
        end

        k_th=k_th+GaussWt_b(1,i)*detJ_b(1,i)*(Bth'*N_th*Bth);
        
        
%-------------- Aerodynamic damping and stiffness matrix---------------------
 
        BA=zeros(ele_dof,5);
        BA_x=zeros(5,ele_dof);
        BA_y=zeros(5,ele_dof);
        
        for j=1:ele_node
            BA(4*ele_node+j,5)=H(j,i);
            BA_x(5,4*ele_node+j)=(dy_deta(1,i)*dH_dxi(j,i)-dy_dxi(1,i)*dH_deta(j,i))/detJ_b(1,i);
            BA_y(5,4*ele_node+j)=(-dx_deta(1,i)*dH_dxi(j,i)+dx_dxi(1,i)*dH_deta(j,i))/detJ_b(1,i);
        end
               
        a_f=a_f+GaussWt_b(1,i)*detJ_b(1,i)*(BA*BA_x*C_phi+BA*BA_y*S_phi);
        a_d=a_d+GaussWt_b(1,i)*detJ_b(1,i)*(BA*BA');

    end

%-------------------------------------------------------------------   
%                         2nd order integration
%------------------------------------------------------------------
%    ------------------- Shear stiffness matrix-------------------

    ks=zeros(ele_dof,ele_dof);     %initialization of ks
 
    dx_dxi=zeros(1,ords);     dx_deta=zeros(1,ords);     dy_dxi=zeros(1,ords);     dy_deta=zeros(1,ords);       %initialize summation terms
    
    detJ_s=zeros(1,ords);            %initialize summation terms
     
    
    for i=1:ords 
        for j=1:ele_node
           dx_dxi(1,i)= dx_dxi(1,i)+X_coord(j)*dHs_dxi(j,i);
           dy_deta(1,i)= dy_deta(1,i)+Y_coord(j)*dHs_deta(j,i);
           dx_deta(1,i)= dx_deta(1,i)+X_coord(j)*dHs_deta(j,i);
           dy_dxi(1,i)= dy_dxi(1,i)+Y_coord(j)*dHs_dxi(j,i);
        end
        
       detJ_s(1,i)=dx_dxi(1,i)*dy_deta(1,i)-dx_deta(1,i)*dy_dxi(1,i);    % Get Jacobian %
       
        Bs=zeros(2,ele_dof);            % initialization of kinematic matrix for bending
        Bs1=zeros(2,ele_dof);
        Bs2=zeros(2,ele_dof);
        
        for j=1:ele_node
            Bs1(1,4*ele_node+j)=(-dHs_dxi(j,i)*dx_deta(1,i)+dHs_deta(j,i)*dx_dxi(1,i))/detJ_s(1,i);
            Bs1(2,4*ele_node+j)=(dHs_dxi(j,i)*dy_deta(1,i)-dHs_deta(j,i)*dy_dxi(1,i))/detJ_s(1,i);
            Bs2(1,ele_node+j)=Hs(j,i);
            Bs2(2,j)=Hs(j,i);       
        end
        Bs=Bs1+Bs2;
            
        ks=ks+GaussWt_s(1,i)*detJ_s(1,i)*(Bs'*As*Bs)*kappa;       
    end

    %_______________summation of kb and ks________
    k_ele_e=zeros(ele_dof,ele_dof);
    k_ele_e=kb+ks;             
   
   
%%%%% Assemble element stiffness matrix to Global stiffness matrix %%%%%%%%   
  
   for i=1:ele_node
       for j=1:ele_node                   
           for p=1:node_dof               
               for q=1:node_dof
                   r=ele_node*(p-1)+i;                          % row index of k_ele matrix
                   s=ele_node*(q-1)+j;                          % column index of k_ele matrix
                   indx1=tot_node*(p-1)+node_indx(ne,i);        % row index of K_global matrix
                   indx2=tot_node*(q-1)+node_indx(ne,j);        % row index of K_global matrix    
                   K_elastic(indx1,indx2)=K_elastic(indx1,indx2)+k_ele_e(r,s);  % elastic stiffness matrix
                   K_thermal(indx1,indx2)=K_thermal(indx1,indx2)+k_th(r,s);     % thermal stiffness matrix
                   A_F(indx1,indx2)=A_F(indx1,indx2)+a_f(r,s);                  % aerodynamic stiffness matrix
                   A_D(indx1,indx2)=A_D(indx1,indx2)+a_d(r,s);                  % aerodynamic damping matrix
               end
           end
       end
   end
   
end