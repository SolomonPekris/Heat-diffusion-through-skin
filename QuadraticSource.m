function [src] = QuadraticSource(f,eID,msh)
%%QUADRATICSOURCE calculates the 3-by-1 source vector in a diffusion reaction equation
%%using quadratic basis functions
%
%   [SRC] = QUADRATICSOURCE(F,EID,MSH) creates a source vector F for a given
%   element ID EID and 1D MSH by using quadratic basis functions, and a
%   gauss scheme with N = 2

gq = CreateGQScheme(2); % Create Gauss scheme
src = 0; % Initialise source vector
J = msh.elem(eID).J; % Set Jacobian

for i = 1:2
    xi = gq.xipts(i); % set xi to gauss point
    psi = [xi*(xi-1)/2; 1-xi^2; xi*(1+xi)/2]; % set psi to quadratic vector
    src = src + gq.gsw(i) * psi; % add gauss quadrature to numerically integrate
end

src = src * f * J; % multiply vector by Jacobian and source