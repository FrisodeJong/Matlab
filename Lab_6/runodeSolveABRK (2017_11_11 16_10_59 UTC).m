clear all; close all
R = 1; P=1; p=5; timeRK=zeros(8,1); timeAB=zeros(8,1);

omega=sqrt((4*pi^2)/R^3);
tRange=[0 5*P];

xexact = @(t) R*[cos(omega.*t); sin(omega.*t); 0]; uexact = @(t) omega*R*[-sin(omega.*t); cos(omega.*t); 0];
solExact = @(t) [uexact(t); xexact(t)];
solVec=solExact(tRange(1));

f = @(t,solVec) twoBodyF(t, solVec); Jac= @twoBodyJac;
h=zeros(8,1);solStruct={};tRK=0;tAB=0;

figerrRK=figure('Name','errorRK','NumberTitle','off'); hold on
figerrAB=figure('Name','errorAB','NumberTitle','off'); hold on

% for l=1:100
%     tic
%     odeSolveRK(f, tRange, solVec, p, 2^(-8)*P);
%     tRK=toc+tRK;
% end
% for l=1:100
%     tic
%     odeSolveAB(f, tRange, solVec, p, 2^(-8)*P);
%     tAB=toc+tAB;
% end
fprintf('Runge-Kutta takes %9.9f seconds to run 100 times \n',tRK)
fprintf('Adams-Bashforth takes %9.9f seconds to run 100 times \n',tAB)
for i=1:8
    h(i)=2^(-i) * P;
    [tArrayRK, solArrayRK] = odeSolveRK(f, tRange, solVec, p, h(i));
    [tArrayAB, solArrayAB] = odeSolveAB(f, tRange, solVec, p, h(i));
    solStruct{i,1}=solArrayRK;
    solStruct{i,2}=solArrayAB;
    for k=1:3
        theta=(k-1)*(1/2);
        [tArrayTheta, solArrayTheta] = odeSolveTheta(f, tRange, solVec, Jac, theta, h(i));
    end
    for k=1:length(tArrayRK)
        errorRK(k,i)=norm(solExact(tArrayRK(k))'-solArrayRK(k,:));
        errorAB(k,i)=norm(solExact(tArrayRK(k))'-solArrayAB(k,:));
    end
    errorEnd(i,1)=errorRK(end,i);
    errorEnd(i,2)=errorAB(end,i);
    
    figure(figerrRK);
    semilogy(tArrayRK, errorRK(:,i));
    Legend{i}=strcat('i=', num2str(i));
    figure(figerrAB);
    semilogy(tArrayAB, errorAB(:,i));
%     errorTheta(i)=norm(solExact-solArrayTheta);
end
figure(figerrRK); legend(Legend); figure(figerrAB); legend(Legend); 

figerrRK=figure('Name','errorEnd','NumberTitle','off'); hold on
plot(h,errorEnd);
legend('ErrorRK at end','ErrorAB at end');
%loglog(erroTheta)