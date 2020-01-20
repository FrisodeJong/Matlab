function [L,U] = LUfactorization(A)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[n1,n2]=size(A);

if n1 ~= n2
    error('Matrix A must be a square matrix')
end

L=zeros(n1,n2);
U=zeros(n1,n2);

for i=1:n1
    for j=1:n2
        L(i,j)=A(i,j);
    end
end

for i=1:n1
    for j=1:n2
        U(i,j)=A(i,j);
    end
end



%outputArg1 = inputArg1;
%outputArg2 = inputArg2;
end

