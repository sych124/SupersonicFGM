function [x] = newmark(M,C,K,F,time,x0,xd0)

% M, C, K are matrices multiplying xddot, xdot and x repectively
% F is column vector of exciations. x0, xd0 are initial x0 and xd vectors

n = size(transpose(time));
t_interv = (time(2)-time(1));

x(:,1) = x0';
xd(:,1) = xd0';

gamma = 1/2; 
beta = 1/4;

A = (1/(beta*t_interv^2))*M+(gamma/(beta*t_interv))*C+K; 
invA = inv(A);

xdd(:,1) = inv(M)*(F(:,1)-C*xd(:,1)-K*x(:,1));

for i = 1:n-2

 B = (F(:,i+1)+ M*((1/(beta*t_interv^2))*x(:,i)+(1/(beta*t_interv))*xd(:,i)+ ...
    (1/(2*beta)-1)*xdd(:,i))+C*((gamma/(beta*t_interv))*x(:,i)+...
    (gamma/beta-1)*xd(:,i)+(gamma/beta-2)*(t_interv/2)*xdd(:,i)));

 x(:,i+1) = invA*B;
 xdd(:,i+1) = (1/(beta*t_interv^2))*(x(:,i+1)-x(:,i))...
    -(1/(beta*t_interv))*xd(:,i)-((1/(2*beta))-1)*xdd(:,i);
 xd(:,i+1) = xd(:,i)+(1-gamma)*t_interv*xdd(:,i)+gamma*t_interv*xdd(:,i+1);

end

x = transpose(x);

end