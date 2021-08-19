%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Thermal_Load.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_thermal=zeros(tot_dof,1);              % thermal load vector for Global system

for ne=1:tot_ele                                       
    f_th=zeros(ele_dof,1);

    dx_dxi=zeros(1,ordb);     dx_deta=zeros(1,ordb);     dy_dxi=zeros(1,ordb);     dy_deta=zeros(1,ordb);       %initialize summation terms
    
    detJ_b=zeros(1,ordb);            %initialize summation terms
    
    for i=1:ordb 
        for j=1:ele_node
           dx_dxi(1,i)= dx_dxi(1,i)+X_coord(j)*dH_dxi(j,i);
           dy_deta(1,i)= dy_deta(1,i)+Y_coord(j)*dH_deta(j,i);
           dx_deta(1,i)= dx_deta(1,i)+X_coord(j)*dH_deta(j,i);
           dy_dxi(1,i)= dy_dxi(1,i)+Y_coord(j)*dH_dxi(j,i);
        end
        
        detJ_b(1,i)=dx_dxi(1,i)*dy_deta(i)-dx_deta(i)*dy_dxi(i);    % Get Jacobian %
       
        B_th1=zeros(ele_dof,3);
        for j=1:ele_node
            B_th1(2*ele_node+j,1)=(dH_dxi(j,i)*dy_deta(1,i)-dH_deta(j,i)*dy_dxi(1,i))/detJ_b(1,i);
            B_th1(2*ele_node+j,3)=(-dH_dxi(j,i)*dx_deta(1,i)+dH_deta(j,i)*dx_dxi(1,i))/detJ_b(1,i);
            B_th1(3*ele_node+j,2)=(-dH_dxi(j,i)*dx_deta(1,i)+dH_deta(j,i)*dx_dxi(1,i))/detJ_b(1,i);
            B_th1(3*ele_node+j,3)=(dH_dxi(j,i)*dy_deta(1,i)-dH_deta(j,i)*dy_dxi(1,i))/detJ_b(1,i);
        end
        
        B_th2=zeros(ele_dof,3);
        for j=1:ele_node
            B_th2(j,1)=(dH_dxi(j,i)*dy_deta(1,i)-dH_deta(j,i)*dy_dxi(1,i))/detJ_b(1,i);
            B_th2(j,3)=(-dH_dxi(j,i)*dx_deta(1,i)+dH_deta(j,i)*dx_dxi(1,i))/detJ_b(1,i);
            B_th2(ele_node+j,2)=(-dH_dxi(j,i)*dx_deta(1,i)+dH_deta(j,i)*dx_dxi(1,i))/detJ_b(1,i);
            B_th2(ele_node+j,3)=(dH_dxi(j,i)*dy_deta(1,i)-dH_deta(j,i)*dy_dxi(1,i))/detJ_b(1,i);
        end
         
        f_th=f_th+GaussWt_b(1,i)*detJ_b(1,i)*(B_th1*N_delta_T+B_th2*M_delta_T);
      
   end
 
  % Assemble f (element) matrix to F_global (Global) matrix
   for i=1:ele_node
       for p=1:node_dof
           r=ele_node*(p-1)+i;                          % row index of f matrix
           indx1=tot_node*(p-1)+node_indx(ne,i);        % row index of F_global matrix
           F_thermal(indx1,1)=F_thermal(indx1,1)+f_th(r,1);
       end
   end
end


