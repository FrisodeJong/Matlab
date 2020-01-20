clear all; close all

mu=10^-3; gamma=pi^2*mu; tEnd=1000; N=20;

p=@(x,t) (pi^2)*mu*sin(pi*x);
uexact=@(x,t) (1+exp(-gamma*t))*sin(pi*x);
u0Func=@(x) uexact(x,0);
g=@(x) 0;

x=linspace(0,1,N+1); xVals=x(2:end-1);
vectheta=[0 0 0 1/2 1/2 1 1];
deltat=[0.625 1.25 1.35 100 200 100 200];

for k=1:7
    dt=deltat(k);
    theta=vectheta(k);
    [tArray, sol, nodes] = heatSolveTheta(p, u0Func , mu, theta, tEnd, N, dt,g);
    error(k)=max(abs(sol(:,end)-uexact(nodes,tEnd)));
    c(k)=mu*deltat(k)/(1/N)^2;
end

table=[vectheta' deltat' c' error'];
array2table(table,'VariableNames', {'theta' 'deltat' 'constant' 'error'})

N=100; dt=1; tEnd=500; theta=0.5;

g=@(x) 293; u0Func=@(x) 293*ones(1,length(x));

sympref('HeavisideAtOrigin', 1);
p=@(x) (heaviside(x-0.4)-heaviside(x-0.6))*2*10^7 /(7.874*10^3 * 4.5*10^2) ;
mu = 8.04*10/(7.874*10^3 * 4.5*10^2);

[tArray, sol, nodes]=heatSolveTheta(p, u0Func , mu, theta, tEnd, N, dt,g);
plot(nodes, sol)

a=find(sol(round(N/2+1),:)>=1811); t=a(1)*dt;
fprintf('bar starts to melt at t=%1.0f',t);