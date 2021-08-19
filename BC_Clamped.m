%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            BC_Clamped.m 
%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:x_node
    bc1(i)=i;                                            % phi_x along lower boundary
    bc2(i)=tot_node+i;                                  % phi_y along lower  boundary
    bc3(i)=2*tot_node+i;                                  % u along lower  boundary
    bc4(i)=3*tot_node+i;                                  % v along lower  boundary
    bc5(i)=4*tot_node+i;                                  % w along lower  boundary
    
    bc6(i)=(y_node-1)*x_node+i;                          % phi_x along upper boundary
    bc7(i)=tot_node+(y_node-1)*x_node+i;                % u along upper boundary
    bc8(i)=2*tot_node+(y_node-1)*x_node+i;                % u along upper boundary
    bc9(i)=3*tot_node+(y_node-1)*x_node+i;                % v along upper boundary
    bc10(i)=4*tot_node+(y_node-1)*x_node+i;                % w along upper boundary
end

for i=1:(y_node-2)
    bc11(i)=x_node*i+1;                                   % phi_x along left side without conner
    bc12(i)=tot_node+x_node*i+1;                         % phi_y along left side without conner
    bc13(i)=2*tot_node+x_node*i+1;                        % u along left side without conner
    bc14(i)=3*tot_node+x_node*i+1;                       % v along left side without conner
    bc15(i)=4*tot_node+x_node*i+1;                       % w along left side without conner
    
    bc16(i)=x_node*(i+1);                                % phi_x along right side without conner
    bc17(i)=1*tot_node+x_node*(i+1);                     % phi_y along right side without conner
    bc18(i)=2*tot_node+x_node*(i+1);                     % u along right side without conner
    bc19(i)=3*tot_node+x_node*(i+1);                     % v along right side without conner
    bc20(i)=4*tot_node+x_node*(i+1);                     % w along right side without conner
    
end


%___________________________________________________________________
bc_dof=sort([bc1 bc2 bc3 bc4 bc5 bc6 bc7 bc8 bc9 bc10 bc11 bc12 bc13 bc14 bc15 bc16 bc17 bc18 bc19 bc20]);
