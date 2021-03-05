function [L2error] = L2(eletype,Ne,theta)
%%L2 calculates L2 error for transient diffusion
%
%   [L2ERROR] = L2(ELETYPE,NE,THETA) calculates the error between exact
%   solutions and numeric solutions. It can consider either "linear" or
%   "quadratic" basis functions using ELETYPE, NE is the number of
%   elements, and THETA is the time stepping method that will be used

%set parameters
x0 = 0;
x1 = 1;
c0 = 0; % dirichlet boundary condition at x = 0
c1 = 1; % dirichlet boundary condition at x = 1
ct0 = 0; % initial concentration state
D = 1; % diffusion coefficient
dt = 0.001; % timestep size
tmax = 1; % end time
msh = OneDimLinearMeshGen(0,1,Ne); % initialise 1D mesh
J = msh.elem(1).J; % set Jacobian

% checks if basis function is linear or quadratic and sets parameters
% accordingly
if eletype == 'linear'
    x = 0:1/Ne:1; % set distance vector
    N = 2 ; % set number of gauss points
elseif eletype == 'quadratic'
    x = 0:1/(Ne*2):1; % set distance vector
    N = 3; % set number of gauss points
end

gq = CreateGQScheme(N); % create gq scheme
c = TransientDiffusion(x0,x1,c0,c1,ct0,D,Ne,dt,tmax,theta,eletype); % calculate numeric soluion across all timesteps
L2error = zeros(1+tmax/dt,1); % initialise L2 error output

for t = 0:dt:tmax % loop through all timesteps
    error = zeros(Ne,1); % initialise error at start of each timestep
    for i = 1:Ne % loop through each elemene
        for ii = 1:N % loop through each gauss point
            xi = gq.xipts(ii); % set gauss value
            if eletype == 'linear' % checks if basis function is linear
                psi = [(1-xi)/2 (1+xi)/2]; % sets psi
                cExact = TransientAnalyticSoln((x(i) * psi(1) + x(i+1) * psi(2)),t); % calculates exact solution
                cNum = c(i,round(t/dt+1))*psi(1)+c(i+1,round(t/dt+1))*psi(2); % calculates numeric solution
            else % use for quadratic basis function
                psi = [xi*(xi-1)/2 (1-xi^2) xi*(1+xi)/2]; %sets psi
                cExact = TransientAnalyticSoln((x(i*2 -1) * psi(1) + x(i*2) * psi(2) + x(i*2 +1)*psi(3)),t); % calculates exact solution
                cNum = c(i*2-1,round(t/dt+1))*psi(1)+c(i*2,round(t/dt+1))*psi(2)+c(i*2+1,round(t/dt+1))*psi(3); % calculates numeric solution
            end
            error(i) = (cExact-cNum)^2 * gq.gsw(ii) + error(i); % sums errors according to gauss scheme
        end
    end
    L2error(round(1+t/dt)) = sqrt(sum(error)*J); % sums all errors for one timestep and multiplies by Jacobian
end
plot(0:dt:tmax,L2error)
% Create ylabel
ylabel({'L2 norm'});

% Create xlabel
xlabel({'Time'});

grid on

set(gca, 'YScale', 'log');