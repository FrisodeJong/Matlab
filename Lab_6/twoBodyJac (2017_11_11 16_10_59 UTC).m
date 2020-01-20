% INPUT
% t         current time (not used)
% solVec    current solution (should be 6x1 array)
% OUTPUT 
% J         Jacobian matrix
function J = twoBodyJac(t, solVec)

B=-4*pi^2/(norm(solVec(1:3))^3)*eye(3)-(12*pi^2)/norm(solVec(1:3))^5*(solVec(1:3)*solVec(1:3)');
J=zeros(6);
J(1:3,4:6)=B;
J(4:6,1:3)=eye(3);