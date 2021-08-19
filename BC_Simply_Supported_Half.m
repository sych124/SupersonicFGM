%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             BC_Simply_Supported.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:x_node
    Bc1(i)=2*tot_node+i;                                  % u along lower  boundary
    Bc2(i)=3*tot_node+i;                                  % v along lower  boundary
    Bc3(i)=4*tot_node+i;                                  % w along lower  boundary
    Bc4(i)=i;                                             % phi_x along lower boundary
    
    Bc5(i)=tot_node+x_node*(y_node-1)+i;                  % phi_y along upper boundary
    Bc6(i)=3*tot_node+x_node*(y_node-1)+i;                % v along upper boundary
end

for i=1:(y_node-2)
    Bc7(i)=2*tot_node+x_node*i+1;                       % u along left side without conner
    Bc8(i)=3*tot_node+x_node*i+1;                       % v along left side without conner
    Bc9(i)=4*tot_node+x_node*i+1;                       % w along left side without conner
    Bc10(i)=2*tot_node+x_node*(i+1);                     % u along right side without conner
    Bc11(i)=3*tot_node+x_node*(i+1);                     % v along right side without conner
    Bc12(i)=4*tot_node+x_node*(i+1);                     % w along right side without conner
    Bc13(i)=1*tot_node+x_node*i+1;                       % phi_y along left side without conner
    Bc14(i)=1*tot_node+x_node*(i+1);                     % phi_y along right side without conner
end

Bc15(1)=1*tot_node+1;                               %phi_y at left lower conner
Bc15(2)=1*tot_node+x_node;                          %phi_y at right lower conner

Bc15(3)=2*tot_node+x_node*(y_node-1)+1;                        %u at left upper conner
Bc15(4)=2*tot_node+x_node*y_node;                              %u at right upper conner

Bc15(5)=3*tot_node+x_node*(y_node-1)+1;                        %v at left upper conner
Bc15(6)=3*tot_node+x_node*y_node;                              %v at right upper conner

Bc15(7)=4*tot_node+x_node*(y_node-1)+1;                        %w at left upper conner
Bc15(8)=4*tot_node+x_node*y_node;                              %w at right upper conner


bc_dof=sort([Bc1 Bc2 Bc3 Bc4 Bc5 Bc6 Bc7 Bc8 Bc9 Bc10 Bc11 Bc12 Bc13 Bc14 Bc15]); 

