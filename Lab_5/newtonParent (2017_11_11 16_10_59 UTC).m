% INPUT
% f         function of rootfinding problem
% df        function returning the Jacobian
% x0        initial guess
% tol       desired tolerance
% maxIt     maximum number of iterations
% OUTPUT    
% root      if succesful, root is an approximation to the
%           root of the nonlinear equation
% flag      flag = 0: tolerance attained, flag = 1: reached maxIt
% iter      the number of iterations
% convHist  convergence history
function [root, flag, convHist] = newtonParent(f, df, x0, tol,...
    maxIt)

convHist=zeros(maxIt,1);
x=x0;
flag=1;
for k=1:maxIt
    xnew=x-(f(x)/df(x));
    convHist(k)=norm(xnew-x);
    if norm(x-xnew)< tol
        flag=0;
        break
    end
    x=xnew;
end
root=x;