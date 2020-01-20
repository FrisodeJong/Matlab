% Solves the 1D heat eqn 
%       du/dt = mu laplace u + p(x,t) on (0,1) x (0, tEnd)
% and initial condition u0 (the Dirichlet boundary 
% conditions are imposed by u0)
% INPUT
% p         right-hand side forcing term (function of x and t)
% u0Func    initial condition function (function of x)
% mu        diffusion coefficient
% theta     parameter for time integration
% tEnd      end time
% N         number of subintervals
% dt        step-size
% OUTPUT
% tArray    array containing the time points
% solArray  array containing the solution at each time level
%           (the ith row equals the solution at time tArray(i))
%           (nrTimeSteps + 1) x (N+1) array
% nodes     (N+1) x 1 array with location of spatial nodes
% g         boundary condition
function [tArray, sol, nodes] = heatSolveTheta(p,...
    u0Func , mu, theta, tEnd, N, dt,g)

x=linspace(0,1,N+1); xVals=x(2:end-1);

h=1/N;
A = makeLaplaceTwo(N);
u0 = u0Func(xVals)'; 
tRange=[0;tEnd];
b=p(xVals)';
b(1) = b(1) + mu/(h^2)*g(0);
b(end) = b(end) + mu/(h^2)*g(1);
f=@(t,u) (-mu/h^2)*A*u + b;
df=@(t,u) (-mu/h^2)*A;

[tArray, solArray] = odeSolveTheta(f, tRange, u0, df, theta, dt);

sol=zeros(N+1,length(solArray(1,:)));
sol(1,:)=g(0)*ones(1,length(solArray(1,:))); 
sol(end,:)=g(1)*ones(1,length(solArray(1,:))); 

sol(2:end-1,:)=solArray;

nodes=x';
end


