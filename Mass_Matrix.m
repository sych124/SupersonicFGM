%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Mass_Matrix.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% You can get Mass matrix from this

M_global=zeros(tot_dof,tot_dof);        % elastic stiffness matrix


i0=I0;
i1=I1;
i2=I2;



ii=zeros(5,5);
ii=[ i2, 0,  i1, 0,  0;
     0,  i2, 0,  i1, 0;
     i1, 0,  i0, 0,  0;
     0,  i1, 0,  i0, 0;
     0,  0,  0,  0,  i0];

for ne=1:tot_ele                                       

    m_ele=zeros(ele_dof,ele_dof);       %initialize element mass matrix
 
    for i=1:ele_node
        X_coord(i)=coord(node_indx(ne,i),1);
        Y_coord(i)=coord(node_indx(ne,i),2);
    end
%-------------------------------------------------------------------  
%                           3rd order integration
%-------------------------------------------------------------------
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
        
        
%    -------------------  element mass matrix-------------------
        Dm=zeros(5,ele_dof);            

       for j=1:ele_node
            Dm(1,j)=H(j,i);
            Dm(2,ele_node+j)=H(j,i);
            Dm(3,2*ele_node+j)=H(j,i);
            Dm(4,3*ele_node+j)=H(j,i);
            Dm(5,4*ele_node+j)=H(j,i);
        end

        m_ele=m_ele+GaussWt_b(1,i)*detJ_b(1,i)*(Dm'*ii*Dm);
        
    end
%%%%% Assemble element stiffness matrix to Global stiffness matrix %%%%%%%%   
  
   for i=1:ele_node
       for j=1:ele_node                   
           for p=1:node_dof               
               for q=1:node_dof
                   r=ele_node*(p-1)+i;                          % row index of k_ele matrix
                   s=ele_node*(q-1)+j;                          % column index of k_ele matrix
                   indx1=tot_node*(p-1)+node_indx(ne,i);        % row index of K_global matrix
                   indx2=tot_node*(q-1)+node_indx(ne,j);        % row index of K_global matrix    
                   M_global(indx1,indx2)=M_global(indx1,indx2)+m_ele(r,s);  % elastic stiffness matrix
                end
           end
       end
   end
   
end