function [root, flag, convHist] = newton(vecF, J, X0, tol,...
    maxIt)
% INPUT
% f         is a string with the name of the function, or a 
%           function handle
% df        is a string with the name of the function, or a 
%           function handle, that can evaluate the Jacobian 
%           matrix of f at a point (or an approximation)
% x0        is the initial guess
% tol       is the desired tolerance
% maxIt     maximum number of iterations
% varargin  parameters for f, df
% OUTPUT    
% root      if succesful, root is an approximation to the
%           root of the nonlinear equation
% flag      flag = 0: tolerance attained, flag = 1: reached maxIt
% iter      the number of iterations
% convHist  convergence history

convHist=zeros(maxIt,1);
[n,n]=size(X0);
X=reshape(X0, n^2, 1);
flag=1;
for k=1:maxIt
    Xnew=X-(inv(J(X))*feval(vecF,X));
    convHist(k)=norm(Xnew-X);
    if norm(X-Xnew)< tol
        flag=0;
        break
    end
    X=Xnew;
end
X=reshape(X,n,n);
root=X;

