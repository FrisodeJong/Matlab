% Attempts to solve A * x = b, with initial guess x0 using 
% an iterative method of the form
%           x^{k+1} = x^{k} + alpha_k P \ r^{k}, r^{k} = b - A x^{k}
% where alpha_k = alpha0 if dynamic = 0
% INPUT
% A         n x n matrix
% b         n x 1 right-hand side
% x0        initial guess
% tol       desired tolerance
% maxIt     maximum number of iterations
% P         preconditioner of A (optional)
% dynamic   0: static, 1: minimise A-norm of error, 2: minimise
%           2-norm of residual
% alpha0    if dynamic = 0, this value is used for alpha_k
% OUTPUT
% x         approximate solution to A * x = b
% flag      if 0 then tolerance is attained
% convHist  relative residual per iteration
function [x, flag, convHist] = iterMethod(A, b, x0, tol, maxIt,...
    P, dynamic, alpha0)

r=b-A*x0;
flag = 1;
x=x0;
if isempty(P)==1;
    P=eye(size(A,1));
end
for k=1:maxIt
    if norm(r) >= tol*norm(b)
    z = P\r;
    if dynamic == 0
    alpha = alpha0;
elseif dynamic == 1
    alpha = (z'*r) / (z'*A*z);
elseif dynamic == 2
    alpha = ((A*z)'*r) / ((A*z)'*A*z);
end
    x = x + alpha*z;
    r = r - alpha*A*z;
    convHist(k)=norm(r);
    else
            flag = 0;
            break
    end
     
end