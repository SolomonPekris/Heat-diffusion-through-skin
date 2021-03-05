function [elemat] = QuadraticDiffusionElem(D,eID,msh)
%%QUADRATICDIFFUSIONELEM calculates the 3-by-3 quadratic local element
%%matrix for the diffusion term of an element in an FEM
%
%   [ELEMAT] = QUADRATICDIFFUSIONELEM(D,EID,MSH) gives the 3-by-3 quadratic
%   local element matrix ELEMAT when given the diffusion coefficient D, the
%   element ID EID, and a 1 dimensional mesh MSH.

gq = CreateGQScheme(3); % Create gauss scheme for N = 3
elemat=0; % Initiate solution matrix
J = msh.elem(eID).J; % Set Jacobian

for n = 1:3 % Loop over Gauss points
 xi = gq.xipts(n); % set xi to Gauss point
 dpsim = [xi-1/2; -2*xi; 1/2+xi]; % Differential vector of psim with respect to xi
 dpsin = [(xi-1/2) (-2*xi) (1/2+xi)]; % Differential vector of psin with respect to xi
 elemat = elemat + gq.gsw(n) * dpsim * dpsin; % sum Gauss quadrants
end

elemat = elemat * D/J; % multiply matrix by coefficients to complete integration
