% Performs integration of the system of ODEs given by
%           d/dt u = f(t, u(t)), u(tRange(1)) = u0
% using the theta-method
% INPUT
% f         the right-hand side function, the output should be a
%           N x 1 array where N is the number of unknowns
% tRange    the time interval of integration
% u0        initial solution at t = tRange(1) (N x 1 array)
% df        a function that evaluates the jacobian of f
% theta     defines the method
% h         the step-size
% OUTPUT
% tArray    array containing the time points
% solArray  array containing the solution at each time level
%           (the ith row equals the solution at time tArray(i))
function [tArray, solArray] = odeSolveThetaTwo(f, tRange, u0, df,...
    theta, h)


u=u0;
tArray=[tRange(1)];

while tArray(end)<tRange(2) 
    new = tArray(end) + h;
    tArray = [tArray new'];  
end

a=size(tArray,2);
n=size(u,1);
I=eye(n);
uN(1,:)=u0;
solArray=zeros(a,n);
solArray(1,:)=u0';

for i=1:a-1
    if theta==0
        u = u + h * f(u,tArray(i));
    else
        F=@(x) u - x + h * ( ...
             theta * f(x, tArray(i + 1)) + (1 - theta) * f(u, tArray(i)) ...
        ); %function uN
     
        J=@(x) theta * h * df(x) - I; %Jacobian for Theta 1/2 and 1
        [u] = newton_pcode(F,J, u, 10^(-8),25); %root= u_n    
    end
    solArray(i+1,:)=u';
end
end
