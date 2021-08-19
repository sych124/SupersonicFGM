%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Static_Nonlinear_Matrix.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% You can get Nonlinear stiffness matrix N1_s & N2_s

N1_s=zeros(tot_dof,tot_dof);        % 1st order Nonlinear stiffness matrix
N2_s=zeros(tot_dof,tot_dof);        % 1st order Nonlinear stiffness matrix

for ne=1:tot_ele                                       
    
    ele_d=zeros(ele_dof,1);
    
    for i=1:ele_node
        X_coord(i)=coord(node_indx(ne,i),1);
        Y_coord(i)=coord(node_indx(ne,i),2);
        
        for j=1:5
            ele_d(ele_node*(j-1)+i,1)=disp_s(tot_node*(j-1)+node_indx(ne,i),1);
        end
    end
%-------------------------------------------------------------------  
%                           3rd order integration
%-------------------------------------------------------------------
    n1_s=zeros(ele_dof,ele_dof);
    n2_s=zeros(ele_dof,ele_dof);

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
        
%    -------------------  Static Nonlinear Matrix-------------------
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
       
        
        aph_vec=A*Bm*ele_d;
        aph1=aph_vec(1,1);
        aph2=aph_vec(2,1);
        aph3=aph_vec(3,1);
        
        aph_m=[ aph1, aph3;
                aph3, aph2];
            

        bth_vec=B*Bm*ele_d;
        bth1=bth_vec(1,1);
        bth2=bth_vec(2,1);
        bth3=bth_vec(3,1);
        
        bth_m=[ bth1, bth3;
                bth3, bth2];
            

        G_s=zeros(3,2);            
        B_g=zeros(4,ele_dof*2);
        B_d=zeros(ele_dof*2,2);
                
        for j=1:ele_node
            B_g(1,4*ele_node+j)=dH_dxi(j,i);
            B_g(2,4*ele_node+j)=dH_deta(j,i);
            B_g(3,9*ele_node+j)=dH_dxi(j,i);
            B_g(4,9*ele_node+j)=dH_deta(j,i);

            B_d(4*ele_node+j,1)=ele_d(4*ele_node+j,1);
            B_d(9*ele_node+j,2)=ele_d(4*ele_node+j,1);
        end
        G_s=(R*B_g*B_d)/detJ_b(1,i);
 
        
        B_theta=zeros(2,ele_dof);
        
        for j=1:ele_node
            B_theta(1,4*ele_node+j)=(dH_dxi(j,i)*dy_deta(1,i)-dH_deta(j,i)*dy_dxi(1,i))/detJ_b(1,i);
            B_theta(2,4*ele_node+j)=(-dH_dxi(j,i)*dx_deta(1,i)+dH_deta(j,i)*dx_dxi(1,i))/detJ_b(1,i);
        end


        n1_s=n1_s+GaussWt_b(1,i)*detJ_b(1,i)*( Bm'*A*G_s*B_theta + B_theta'*G_s'*A*Bm  + Bb'*B*G_s*B_theta + B_theta'*G_s'*B*Bb + B_theta'*( aph_m + bth_m )*B_theta );
        n2_s=n2_s+GaussWt_b(1,i)*detJ_b(1,i)*(3/2*B_theta'*G_s'*A*G_s*B_theta);
        
    end    
   
%%%%% Assemble element stiffness matrix to Global stiffness matrix %%%%%%%%   
  
   for i=1:ele_node
       for j=1:ele_node
           for p=1:node_dof
               for q=1:node_dof
                   r=ele_node*(p-1)+i;                          % row index of k_ele matrix
                   s=ele_node*(q-1)+j;                          % column index of k_ele matrix
                   indx1=tot_node*(p-1)+node_indx(ne,i);        % row index of K_global matrix
                   indx2=tot_node*(q-1)+node_indx(ne,j);        % column index of K_global matrix    
                   N1_s(indx1,indx2)=N1_s(indx1,indx2)+n1_s(r,s);  % elastic stiffness matrix
                   N2_s(indx1,indx2)=N2_s(indx1,indx2)+n2_s(r,s);     % thermal stiffness matrix
               end
           end
       end
   end
   
end
