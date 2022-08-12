function df = gradobj(theta2)
% This function computes the gradient of the objective function

% Written by Aviv Nevo, May 1998.

global invA IV 

load gmmresid 
load mvalold
temp = jacob(mvalold,theta2)';
df = 2*temp*IV*invA*IV'*gmmresid