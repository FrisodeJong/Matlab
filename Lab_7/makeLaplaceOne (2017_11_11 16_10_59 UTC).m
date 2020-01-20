function L = makeLaplaceOne(N)
% INPUT
% N         number of subintervals
% OUTPUT
% L         discrete Laplace operator (3-point stencil)

a=zeros(N+1,N+1);
for k=2:N
    a(k,k-1)=-1;
    a(k,k)=2;
    a(k,k+1)=-1;
end
a(1,1)=1;
a(N+1,N+1)=1;
L=a;