%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            BC_Clamped.m 
%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:x_node
    Bc1(i)=i;                                            % phi_x along lower boundary
    Bc2(i)=tot_node+i;                                  % phi_y along lower  boundary
    Bc3(i)=2*tot_node+i;                                  % u along lower  boundary
    Bc4(i)=3*tot_node+i;                                  % v along lower  boundary
    Bc5(i)=4*tot_node+i;                                  % w along lower  boundary
    
    Bc6(i)=tot_node+(y_node-1)*x_node+i;                % phi_y along upper boundary
    Bc7(i)=3*tot_node+(y_node-1)*x_node+i;                % v along upper boundary
end

for i=1:(y_node-2)
    Bc8(i)=x_node*i+1;                                   % phi_x along left side without conner
    Bc9(i)=tot_node+x_node*i+1;                         % phi_y along left side without conner
    Bc10(i)=2*tot_node+x_node*i+1;                        % u along left side without conner
    Bc11(i)=3*tot_node+x_node*i+1;                       % v along left side without conner
    Bc12(i)=4*tot_node+x_node*i+1;                       % w along left side without conner
    
    Bc13(i)=x_node*(i+1);                                % phi_x along right side without conner
    Bc14(i)=1*tot_node+x_node*(i+1);                     % phi_y along right side without conner
    Bc15(i)=2*tot_node+x_node*(i+1);                     % u along right side without conner
    Bc16(i)=3*tot_node+x_node*(i+1);                     % v along right side without conner
    Bc17(i)=4*tot_node+x_node*(i+1);                     % w along right side without conner
    
end

Bc18(1)=x_node*(y_node-1)+1;                        %phi_x at left upper conner
Bc18(2)=x_node*y_node;                              %phi_x at right upper conner

Bc18(3)=2*tot_node+x_node*(y_node-1)+1;                        %u at left upper conner
Bc18(4)=2*tot_node+x_node*y_node;                              %u at right upper conner

Bc18(5)=4*tot_node+x_node*(y_node-1)+1;                        %w at left upper conner
Bc18(6)=4*tot_node+x_node*y_node;                              %w at right upper conner

%___________________________________________________________________
bc_dof=sort([Bc1 Bc2 Bc3 Bc4 Bc5 Bc6 Bc7 Bc8 Bc9 Bc10 Bc11 Bc12 Bc13 Bc14 Bc15 Bc16 Bc17 Bc18]); 
