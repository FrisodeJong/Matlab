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
function [tArray, solArray] = odeSolveTheta(f, tRange, u0, df,...
    theta, h)

tArray=[tRange(1):h:tRange(2)];
[m,n]=size(tArray);
solArray(1,:)=u0;
if theta==0
    for t=1:n-1
        solArray(t+1,:)=solArray(t)+h*f(tArray(t),solArray(t,:));
    end
else
    for t=1:n-1
        F = @(u) solArray(t,:) + h*(theta*f(tArray(t+1),u)+(1-theta)*f(tArray(t),solArray(t,:))) - u;
        J= @(u) theta*h*df(t,u)-eye(length(u0));
        solArray(t+1,:)=newton_pcode(F, J, u0,10^-9,25);
    end
end
tArray=tArray';
solArray=solArray';