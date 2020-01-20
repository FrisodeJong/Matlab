tRange=[0;10];
u0=1; tzero=0;thalf=0;tone=0;
lambda=-100;
f = @(t,u) lambda*(u - sin(t)) + cos(t); 
df = lambda;

solutions = {};
for i=1:7
    h=2*10^-2 + (i-4)*10^-4;
    for k=1:3
        theta=(k-1)*(1/2);
        [tArraytemp, solArraytemp] = odeSolveTheta(f, tRange, u0, df, theta, h);
        sol100{k,i}=solArraytemp;
        t100{k, i} = tArraytemp;
    end
end
figure('Name','lambda=100 theta=0','NumberTitle','off')
%for i>4 no convergence
for i=1:7
    plot(t100{1,i},sol100{1,i})
    hold on
end
h=2*10^-2 + (3-4)*10^-4;
fprintf('For i=3 abs(h*lambda+1) equals %9.9f \n',abs(h*lambda+1))
h=2*10^-2 + (4-4)*10^-4;
fprintf('For i=4 abs(h*lambda+1) equals %9.9f \n',abs(h*lambda+1))
legend('1','2','3','4','5','6','7')

figure('Name','lambda=100 theta=0.5','NumberTitle','off')
for i=1:7
    plot(t100{2,i},sol100{2,i})
    hold on
end
legend('1','2','3','4','5','6','7')
figure('Name','lambda=100 theta=1','NumberTitle','off')
for i=1:7
    plot(t100{3,i},sol100{3,i})
    hold on
end
legend('1','2','3','4','5','6','7')
hold off

lambda=-1;
f = @(t,u) lambda*(u - sin(t)) + cos(t); 
df = lambda;
uSolexact = @(t) exp(-1*t) + sin(t);

for i=1:8
    h(i)=2^-(i-1);
     for k=1:3
        theta=(k-1)*(1/2);
        [tArraytemp, solArraytemp] = odeSolveTheta(f, tRange, u0, df, theta, h(i));
        sol1{k,i}=solArraytemp;
        t1{k, i} = tArraytemp;
        error(i,k) = norm(solArraytemp(end)-uSolexact(tRange(end)));
    end
end
for l=1:100
    tic
    odeSolveTheta(f, tRange, u0, df, 0, h(i));
    tzero=toc+tzero;
end
for l=1:100
    tic
    odeSolveTheta(f, tRange, u0, df, 0.5, h(i));
    thalf=toc+thalf;
end
for l=1:100
    tic
    odeSolveTheta(f, tRange, u0, df, 1, h(i));
    tone=toc+tone;
end
fprintf('Forward Euler takes %9.9f seconds to run 100 times \n',tzero)
fprintf('Crank-Nicholson takes %9.9f seconds to run 100 times \n',thalf)
fprintf('Backward Euler takes %9.9f seconds to run 100 times \n',tone)
figure('Name','error','NumberTitle','off')
loglog(h,error(:,1))
hold on
loglog(h,error(:,2))
loglog(h,error(:,3))
legend('errortheta0','errortheta.5','errortheta1')
hold off
