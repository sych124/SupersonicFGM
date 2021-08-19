%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             BC_Simply_Supported.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:x_node
    bc1(i)=2*tot_node+i;                                  % u along lower  boundary
    bc2(i)=3*tot_node+i;                                  % v along lower  boundary
    bc3(i)=4*tot_node+i;                                  % w along lower  boundary
    bc4(i)=2*tot_node+(y_node-1)*x_node+i;                % u along upper boundary
    bc5(i)=3*tot_node+(y_node-1)*x_node+i;                % v along upper boundary
    bc6(i)=4*tot_node+(y_node-1)*x_node+i;                % w along upper boundary
    bc7(i)=i;                                            % phi_x along lower boundary
    bc8(i)=(y_node-1)*x_node+i;                          % phi_x along upper boundary
end

for i=1:(y_node-2)
    bc9(i)=2*tot_node+x_node*i+1;                       % u along left side without conner
    bc10(i)=3*tot_node+x_node*i+1;                       % v along left side without conner
    bc11(i)=4*tot_node+x_node*i+1;                       % w along left side without conner
    bc12(i)=2*tot_node+x_node*(i+1);                     % u along right side without conner
    bc13(i)=3*tot_node+x_node*(i+1);                     % v along right side without conner
    bc14(i)=4*tot_node+x_node*(i+1);                     % w along right side without conner
    bc15(i)=1*tot_node+x_node*i+1;                       % phi_y along left side without conner
    bc16(i)=1*tot_node+x_node*(i+1);                     % phi_y along right side without conner
end

bc17(1)=1*tot_node+1;                               %phi_y at left lower conner
bc17(2)=1*tot_node+x_node;                          %phi_y at right lower conner
bc17(3)=1*tot_node+x_node*(y_node-1)+1;             %phi_y at left upper conner
bc17(4)=1*tot_node+x_node*y_node;                   %phi_y at left upper conner



bc_dof=sort([bc1 bc2 bc3 bc4 bc5 bc6 bc7 bc8 bc9 bc10 bc11 bc12 bc13 bc14 bc15 bc16 bc17]);         

