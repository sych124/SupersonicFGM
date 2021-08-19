function [SCF] = ShrCor(n, E_z, nu,h0,z)
%property


Q11=E_z./(1-nu.^2);
Q12=nu.*E_z./(1-nu.^2);
Q16=0;
Q55=E_z./(2*(1+nu));


A11=0; A12=0; A16=0; A55=0;
D11=0; D12=0; D16=0; D55=0;

dz=z(2)-z(1);

% %  A B D
for i=1:n+1;
    
    A11=A11+Q11(i)*dz;
    A12=A12+Q12(i)*dz;
    A55=A55+Q55(i)*dz;
    
    
    D11=D11+((z(i)-h0)^2)*Q11(i)*dz;
    D12=D12+((z(i)-h0)^2)*Q12(i)*dz;
    D55=D55+((z(i)-h0)^2)*Q55(i)*dz; 
end

A22=A11;
D22=D11;

A=[A11 A12 A16; A12 A22 A16; A16 A16 A55];
B = diag([0 0 0]);
D=[D11 D12 D16; D12 D22 D16; D16 D16 D55];

stiffness=[[A] [B];[B] [D]];
Stiffness=inv(stiffness);

a11=Stiffness(1,1);
b11=Stiffness(1,4);
b12=Stiffness(1,5);
b16=0;
b55=Stiffness(3,6);

d11=Stiffness(4,4);
d12=Stiffness(4,5);
d16=0;
d55=Stiffness(6,6);


SCF = 0;
K=0;
KS=0;

for i=1:n+1
    for j=1:i
        
        KS=KS+((z(j)-h0)*Q11(j)*d11+(z(j)-h0)*Q12(j)*d12)*dz;
    end
    K=K+(KS^2)/Q55(i)*dz;
    KS=0;
end

SCF=(K^(-1))/A55;

end
