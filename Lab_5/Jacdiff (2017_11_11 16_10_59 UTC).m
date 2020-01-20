n=3;
X=eye(n);
vecX=reshape(X, n^2, 1);
B=[2 1 0; 0 2 1; 0 0 2];
J=matJac(vecX, 'exact',1,B);
for l=0:14
    h=10^-l;
    appJ=matJac(vecX, 'approx',h,B);
    A=max(max(appJ-J));
    C=max(max(-(appJ-J)));
    error(l+1)=max(A,C);
end
[m,k]=min(error);
fprintf('minimum error is %s, for h=%d \n', m, 10^-k)

figure
loglog(error)