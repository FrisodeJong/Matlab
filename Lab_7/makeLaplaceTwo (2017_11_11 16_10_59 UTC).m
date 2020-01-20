% INPUT
% N         number of subintervals
% OUTPUT
% L         discrete Laplace operator (3-point stencil)
function L = makeLaplaceTwo(N)
L=2*eye(N-1)-diag(ones(N-2,1),1)-diag(ones(N-2,1),-1);
 
end