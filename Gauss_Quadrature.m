%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Gauss_Quadrature.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%-----------for bending (3 by 3)
Ptb=3;                                              %integration Point
GPtb=[-sqrt(3/5), 0.0, sqrt(3/5)];%[-.906179846,-.53846931,0,.53846931,.906179846];
GWtb=[5/9,8/9,5/9];%[.2369269,.4786287,.5688889,.4786287,.2369269];
for i=1:Ptb
    for j=1:Ptb
        GaussPt_b(Ptb*(j-1)+i,1) = GPtb(1,i);
        GaussPt_b(Ptb*(j-1)+i,2) = GPtb(1,j);
        GaussWt_b(Ptb*(j-1)+i)=GWtb(1,i)*GWtb(1,j);
    end
end

%-----------for bending (2 by 2)
Pts=2;
GPts = [-0.577350269,0.577350269];
GWts = [1.000, 1.000];
for i=1:Pts
    for j=1:Pts
        GaussPt_s(Pts*(j-1)+i,1) = GPts(1,i);
        GaussPt_s(Pts*(j-1)+i,2) = GPts(1,j);
        GaussWt_s(Pts*(j-1)+i)=GWts(1,i)*GWts(1,j);
    end
end

ordb = Ptb*Ptb ;         % Gauss integration order of Bending Term for  w.r.t x
ords = Pts*Pts ;         % Gauss integration order of shear Term for  w.r.t y
