% Solves the 1D Poisson problem 
%       -laplace u = f on (0,1)
% with Dirichlet bdy conditions given by g
% INPUT
% f         right-hand side function
% g         Dirichlet boundary condition function
% N         number of subintervals
% precon    'none', 'jacobi' or 'gs'
% tol       tolerance for iterative solver
% OUTPUT
% sol       (N+1) x 1 solution array
% nodes     (N+1) x 1 array with location of spatial nodes
function [sol, nodes] = poissonSolveFD(f, g, N, precon, tol, h)

u(1)=g(0);
u(end)=g(1);
for k=1:N+1
    nodesArray(k)=(k-1)*h;
end
for i=2:N
    fArray(i)=f(nodesArray(i));
end
fArray(1)=f(nodesArray(1))+u(1)/h^2;
fArray(N+1)=f(nodesArray(N+1))+u(end)/h^2;
b=h^2 *fArray';
A = makeLaplaceOne(N);
% if strcmpi(precon,'none')==1
%     P = eye(N+1);
% elseif strcmpi(precon,'jacobi')==1
%     P=diag(diag(A));
% elseif strcmpi(precon,'gs')==1
%     P=triu(A);
% end
% [x, flag,R]=iterMethod(A, b, zeros(N+1,1), tol, 25, P, 2, 0);
x=A\b;
sol=x';
nodes=nodesArray;