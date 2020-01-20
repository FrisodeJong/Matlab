% INPUT
% t         current time (not used)
% solVec    current solution (should be 6x1 array)
% OUTPUT 
% dSolVec   right-hand side
function dSolVec = twoBodyF(t, solVec)

dSolVec(1:3)=(-4*pi^2/(norm(solVec(4:6))^3))*solVec(4:6);
dSolVec(4:6)=solVec(1:3);
dSolVec=dSolVec';
end