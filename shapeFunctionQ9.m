% .............................................................
function [shape,naturalDerivatives]=shapeFunctionQ9(xi,eta)
% shape function and derivatives for Q9 elements
% shape : Shape functions
% naturalDerivatives: derivatives w.r.t. xi
% xi: natural coordinates (-1 ... +1)
% eta: natural coordinates (-1 ... +1)
shape=1/4*[xi*eta*(xi-1)*(eta-1);
        xi*eta*(xi+1)*(eta-1);
        xi*eta*(xi+1)*(eta+1);
        xi*eta*(xi-1)*(eta+1);
        -2*eta*(xi+1)*(xi-1)*(eta-1);
        -2*xi*(xi+1)*(eta+1)*(eta-1);
        -2*eta*(xi+1)*(xi-1)*(eta+1);
        -2*xi*(xi-1)*(eta+1)*(eta-1);
        4*(xi+1)*(xi-1)*(eta+1)*(eta-1)];
naturalDerivatives=...
        1/4*[eta*(2*xi-1)*(eta-1),xi*(xi-1)*(2*eta-1);
        eta*(2*xi+1)*(eta-1),xi*(xi+1)*(2*eta-1);
        eta*(2*xi+1)*(eta+1),xi*(xi+1)*(2*eta+1);
        eta*(2*xi-1)*(eta+1),xi*(xi-1)*(2*eta+1);
        -4*xi*eta*(eta-1), -2*(xi+1)*(xi-1)*(2*eta-1);
        -2*(2*xi+1)*(eta+1)*(eta-1),-4*xi*eta*(xi+1);
        -4*xi*eta*(eta+1), -2*(xi+1)*(xi-1)*(2*eta+1);
        -2*(2*xi-1)*(eta+1)*(eta-1),-4*xi*eta*(xi-1);
        8*xi*(eta^2-1), 8*eta*(xi^2-1)];
end % end function shapeFunctionQ9