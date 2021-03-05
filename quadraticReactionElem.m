function [elemat] = quadraticReactionElem(lamda,eID,msh)
%%QUADRATICREACTIONELEM calculates the local 3-by-3 quadratic element matrix for the linear
%%reaction operator, for any element in the finite element mesh for a given
%%reaction coefficient lamda, for an element eID, in a linear 1D msh using
%%a gauss scheme. For lamda = 1, this is also equal to the local element
%%Mass matrix

gq = CreateGQScheme(3); % Generate gauss scheme
elemat=0; % Initiate solution matrix
J = msh.elem(eID).J; % Set Jacobian

for n=1:3 % Loops through each of the gauss scheme points
 xi = gq.xipts(n); % sets xi to the value of the point
 psim = [xi*(xi-1)/2; 1-xi^2; xi*(1+xi)/2]; % create quadratic psi m vector
 psin = [xi*(xi-1)/2 (1-xi^2) xi*(1+xi)/2]; % create quadratic psi n vector
 elemat = elemat + gq.gsw(n) * psim * psin; % sum gauss scheme values to numerically integrate
end
elemat = elemat * lamda * J; % multiply element by lamda and Jacobian to complete numerical integration

