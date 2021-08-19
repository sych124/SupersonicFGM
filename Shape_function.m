%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Shape_function.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculation of the value of the shape function at Gauss points for bending

for i=1:ordb
    H(1,i)=1.0/4.0*(GaussPt_b(i,1)^2-GaussPt_b(i,1))*(GaussPt_b(i,2)^2-GaussPt_b(i,2));             % the Value of H1 at ith Gauss point
    H(2,i)=1.0/4.0*(GaussPt_b(i,1)^2+GaussPt_b(i,1))*(GaussPt_b(i,2)^2-GaussPt_b(i,2));             % the Value of H1 at ith Gauss point
    H(3,i)=1.0/4.0*(GaussPt_b(i,1)^2+GaussPt_b(i,1))*(GaussPt_b(i,2)^2+GaussPt_b(i,2));             % the Value of H1 at ith Gauss point
    H(4,i)=1.0/4.0*(GaussPt_b(i,1)^2-GaussPt_b(i,1))*(GaussPt_b(i,2)^2+GaussPt_b(i,2));             % the Value of H1 at ith Gauss point
    H(5,i)=1/2*(1-GaussPt_b(i,1)^2)*(GaussPt_b(i,2)^2-GaussPt_b(i,2));                              % the Value of H5 at ith Gauss point
    H(6,i)=1/2*(GaussPt_b(i,1)^2+GaussPt_b(i,1))*(1-GaussPt_b(i,2)^2);                              % the Value of H6 at ith Gauss point    
    H(7,i)=1/2*(1-GaussPt_b(i,1)^2)*(GaussPt_b(i,2)^2+GaussPt_b(i,2));                              % the Value of H7 at ith Gauss point
    H(8,i)=1/2*(GaussPt_b(i,1)^2-GaussPt_b(i,1))*(1-GaussPt_b(i,2)^2);                              % the Value of H8 at ith Gauss point
    H(9,i)=(1-GaussPt_b(i,1)^2)*(1-GaussPt_b(i,2)^2);                                               % the Value of H9 at ith Gauss point
end

for i=1:ordb
    dH_dxi(1,i)=1.0/4.0*(2*GaussPt_b(i,1)-1)*(GaussPt_b(i,2)^2-GaussPt_b(i,2));             % the Value of dH1_dxi at ith Gauss point
    dH_dxi(2,i)=1.0/4.0*(2*GaussPt_b(i,1)+1)*(GaussPt_b(i,2)^2-GaussPt_b(i,2));             % the Value of dH2_dxi at ith Gauss point
    dH_dxi(3,i)=1.0/4.0*(2*GaussPt_b(i,1)+1)*(GaussPt_b(i,2)^2+GaussPt_b(i,2));             % the Value of dH3_dxi at ith Gauss point
    dH_dxi(4,i)=1.0/4.0*(2*GaussPt_b(i,1)-1)*(GaussPt_b(i,2)^2+GaussPt_b(i,2));             % the Value of dH4_dxi at ith Gauss point
    dH_dxi(5,i)=1/2*(-2*GaussPt_b(i,1))*(GaussPt_b(i,2)^2-GaussPt_b(i,2));                  % the Value of dH5_dxi at ith Gauss point
    dH_dxi(6,i)=1/2*(2*GaussPt_b(i,1)+1)*(1-GaussPt_b(i,2)^2);                              % the Value of dH6_dxi at ith Gauss point
    dH_dxi(7,i)=1/2*(-2*GaussPt_b(i,1))*(GaussPt_b(i,2)^2+GaussPt_b(i,2));                  % the Value of dH7_dxi at ith Gauss point
    dH_dxi(8,i)=1/2*(2*GaussPt_b(i,1)-1)*(1-GaussPt_b(i,2)^2);                              % the Value of dH8_dxi at ith Gauss point
    dH_dxi(9,i)=-2*GaussPt_b(i,1)*(1-GaussPt_b(i,2)^2);                                     % the Value of dH9_dxi at ith Gauss point
end

for i=1:ordb
    dH_deta(1,i)=1.0/4.0*(GaussPt_b(i,1)^2-GaussPt_b(i,1))*(2*GaussPt_b(i,2)-1);                     % the Value of dH1_deta at ith Gauss point
    dH_deta(2,i)=1.0/4.0*(GaussPt_b(i,1)^2+GaussPt_b(i,1))*(2*GaussPt_b(i,2)-1);                     % the Value of dH2_deta at ith Gauss point
    dH_deta(3,i)=1.0/4.0*(GaussPt_b(i,1)^2+GaussPt_b(i,1))*(2*GaussPt_b(i,2)+1);                     % the Value of dH3_deta at ith Gauss point
    dH_deta(4,i)=1.0/4.0*(GaussPt_b(i,1)^2-GaussPt_b(i,1))*(2*GaussPt_b(i,2)+1);                     % the Value of dH4_deta at ith Gauss point
    dH_deta(5,i)=1/2*(1-GaussPt_b(i,1)^2)*(2*GaussPt_b(i,2)-1);                                      % the Value of dH5_deta at ith Gauss point
    dH_deta(6,i)=1/2*(GaussPt_b(i,1)^2+GaussPt_b(i,1))*(-2*GaussPt_b(i,2));                          % the Value of dH6_deta at ith Gauss point
    dH_deta(7,i)=1/2*(1-GaussPt_b(i,1)^2)*(2*GaussPt_b(i,2)+1);                                      % the Value of dH7_deta at ith Gauss point
    dH_deta(8,i)=1/2*(GaussPt_b(i,1)^2-GaussPt_b(i,1))*(-2*GaussPt_b(i,2));                          % the Value of dH8_deta at ith Gauss point
    dH_deta(9,i)=(1-GaussPt_b(i,1)^2)*(-2*GaussPt_b(i,2));                                           % the Value of dH9_deta at ith Gauss point
end

% Calculation of the value of the shape function at Gauss points For Shear
for i=1:ords
    Hs(1,i)=1.0/4.0*(GaussPt_s(i,1)^2-GaussPt_s(i,1))*(GaussPt_s(i,2)^2-GaussPt_s(i,2));             % the Value of Hs1 at ith Gauss point
    Hs(2,i)=1.0/4.0*(GaussPt_s(i,1)^2+GaussPt_s(i,1))*(GaussPt_s(i,2)^2-GaussPt_s(i,2));             % the Value of Hs1 at ith Gauss point
    Hs(3,i)=1.0/4.0*(GaussPt_s(i,1)^2+GaussPt_s(i,1))*(GaussPt_s(i,2)^2+GaussPt_s(i,2));             % the Value of Hs1 at ith Gauss point
    Hs(4,i)=1.0/4.0*(GaussPt_s(i,1)^2-GaussPt_s(i,1))*(GaussPt_s(i,2)^2+GaussPt_s(i,2));             % the Value of Hs1 at ith Gauss point
    Hs(5,i)=1/2*(1-GaussPt_s(i,1)^2)*(GaussPt_s(i,2)^2-GaussPt_s(i,2));                              % the Value of Hs5 at ith Gauss point
    Hs(6,i)=1/2*(GaussPt_s(i,1)^2+GaussPt_s(i,1))*(1-GaussPt_s(i,2)^2);                              % the Value of Hs6 at ith Gauss point    
    Hs(7,i)=1/2*(1-GaussPt_s(i,1)^2)*(GaussPt_s(i,2)^2+GaussPt_s(i,2));                              % the Value of Hs7 at ith Gauss point
    Hs(8,i)=1/2*(GaussPt_s(i,1)^2-GaussPt_s(i,1))*(1-GaussPt_s(i,2)^2);                              % the Value of Hs8 at ith Gauss point
    Hs(9,i)=(1-GaussPt_s(i,1)^2)*(1-GaussPt_s(i,2)^2);                                               % the Value of Hs9 at ith Gauss point
end


for i=1:ords
    dHs_dxi(1,i)=1.0/4.0*(2*GaussPt_s(i,1)-1)*(GaussPt_s(i,2)^2-GaussPt_s(i,2));             % the Value of dHs1_dxi at ith Gauss point
    dHs_dxi(2,i)=1.0/4.0*(2*GaussPt_s(i,1)+1)*(GaussPt_s(i,2)^2-GaussPt_s(i,2));             % the Value of dHs2_dxi at ith Gauss point
    dHs_dxi(3,i)=1.0/4.0*(2*GaussPt_s(i,1)+1)*(GaussPt_s(i,2)^2+GaussPt_s(i,2));             % the Value of dHs3_dxi at ith Gauss point
    dHs_dxi(4,i)=1.0/4.0*(2*GaussPt_s(i,1)-1)*(GaussPt_s(i,2)^2+GaussPt_s(i,2));             % the Value of dHs4_dxi at ith Gauss point
    dHs_dxi(5,i)=1/2*(-2*GaussPt_s(i,1))*(GaussPt_s(i,2)^2-GaussPt_s(i,2));                  % the Value of dHs5_dxi at ith Gauss point
    dHs_dxi(6,i)=1/2*(2*GaussPt_s(i,1)+1)*(1-GaussPt_s(i,2)^2);                              % the Value of dHs6_dxi at ith Gauss point
    dHs_dxi(7,i)=1/2*(-2*GaussPt_s(i,1))*(GaussPt_s(i,2)^2+GaussPt_s(i,2));                  % the Value of dHs7_dxi at ith Gauss point
    dHs_dxi(8,i)=1/2*(2*GaussPt_s(i,1)-1)*(1-GaussPt_s(i,2)^2);                              % the Value of dHs8_dxi at ith Gauss point
    dHs_dxi(9,i)=-2*GaussPt_s(i,1)*(1-GaussPt_s(i,2)^2);                                     % the Value of dHs9_dxi at ith Gauss point
end


for i=1:ords
    dHs_deta(1,i)=1.0/4.0*(GaussPt_s(i,1)^2-GaussPt_s(i,1))*(2*GaussPt_s(i,2)-1);                         % the Value of dHs1_deta at ith Gauss point
    dHs_deta(2,i)=1.0/4.0*(GaussPt_s(i,1)^2+GaussPt_s(i,1))*(2*GaussPt_s(i,2)-1);                         % the Value of dHs2_deta at ith Gauss point
    dHs_deta(3,i)=1.0/4.0*(GaussPt_s(i,1)^2+GaussPt_s(i,1))*(2*GaussPt_s(i,2)+1);                         % the Value of dHs3_deta at ith Gauss point
    dHs_deta(4,i)=1.0/4.0*(GaussPt_s(i,1)^2-GaussPt_s(i,1))*(2*GaussPt_s(i,2)+1);                         % the Value of dHs4_deta at ith Gauss point
    dHs_deta(5,i)=1/2*(1-GaussPt_s(i,1)^2)*(2*GaussPt_s(i,2)-1);                                          % the Value of dHs5_deta at ith Gauss point
    dHs_deta(6,i)=1/2*(GaussPt_s(i,1)^2+GaussPt_s(i,1))*(-2*GaussPt_s(i,2));                              % the Value of dHs6_deta at ith Gauss point
    dHs_deta(7,i)=1/2*(1-GaussPt_s(i,1)^2)*(2*GaussPt_s(i,2)+1);                                          % the Value of dHs7_deta at ith Gauss point
    dHs_deta(8,i)=1/2*(GaussPt_s(i,1)^2-GaussPt_s(i,1))*(-2*GaussPt_s(i,2));                              % the Value of dHs8_deta at ith Gauss point
    dHs_deta(9,i)=(1-GaussPt_s(i,1)^2)*(-2*GaussPt_s(i,2));                                               % the Value of dHs9_deta at ith Gauss point                                           % the Value of dH4_deta at ith Gauss point
end