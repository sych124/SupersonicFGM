function [P1, P2, result] = MatProp(Pc,Pm,T,k,dmlV,temper_en)
    Vc = (dmlV+1/2).^k; Vm = 1-Vc; % Volume Fraction Matrix
    P2 = Pm(2) + Pm(2) * temper_en*(Pm(1)/T + Pm(3)*T + Pm(4)*T^2 + Pm(5)*T^3);
    P1 = Pc(2) + Pc(2) * temper_en*(Pc(1)/T + Pc(3)*T + Pc(4)*T^2 + Pc(5)*T^3);
    pre_result = P2*Vm + P1*Vc; % material properties
    result = pre_result;
end