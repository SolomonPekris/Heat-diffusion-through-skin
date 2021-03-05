function [c] = TransientDiffusion(x0,x1,c0,c1,ct0,D,Ne,dt,tmax,theta,eletype)
%%TRANSIENTDIFFUSION solves a 1D transient diffusion equation with
%%Dirichlet boundary conditions
%
%   [C] = TRANSIENTDIFFUSION(X0,X1,C0,C1,NE,DT,TMAX,THETA,ELETYPE) solves a
%   transient diffusion equation numerically in 1D. This is between the
%   points X0 and X1, using Dirichlet boundary conditions and these point,
%   C0 and C1 respectively, and ct0 is the initial condition. D is the diffusion coefficient,
%   NE is the number of elements in the mesn, DT is
%   the timestep, and TMAX is the time inteval. The solution can be
%   calculated at each timestep using Forward Euler method, Crank-Nicolson,
%   or Backward Euler, by setting THETA to 0,1/2, or 1 respectively. The
%   basis equation type can be set to "linear" or "quadratic" with ELETYPE.

%set constants
msh = OneDimLinearMeshGen(x0,x1,Ne); % generate 1D mesh

if eletype == "linear" % uses linear basis function
    c(:,1) = zeros(Ne+1,1); % initialise concentration vector
    c(:,1) = ct0; % set initial condition
    x = x0:(x1-x0)/Ne:x1; % Distance vector
    m = zeros(Ne+1,Ne+1); % initialise global mass matrix
    k = zeros(Ne+1,Ne+1); % initialise global stiffness diffusion matrix
    
    for i = 1:Ne % loop through all elements calculating local mass and diffusion matrices  
        % adds local diffusion matrices to global stifness matrix. No reaction
        % or source term term so only diffusion is considered
        k(i:i+1,i:i+1)= k(i:i+1,i:i+1)+ DiffusionElemMatrix(1,i,msh);
        % local mass matrix calculated using the reaction solver with lamda = 1
        m(i:i+1,i:i+1)= m(i:i+1,i:i+1)+ ReactionElemMatrix(1,i,msh);
    end
    
elseif eletype == "quadratic" % usese quadratic basis function
    c(:,1) = zeros(2*Ne+1,1); % initialise concentration vector
    c(:,1) = ct0; % set initial condition    
    x = x0:(x1-x0)/(2*Ne):x1; % Distance vector
    m = zeros(2*Ne+1,2*Ne+1); % initialise global mass matrix
    k = zeros(2*Ne+1,2*Ne+1); % initialise global stiffness diffusion matrix
    
    for i = 1:Ne % loops through each element
    % adds local diffusion matrices to global stifness matrix. No reaction
    % or source term term so only diffusion is considered
    m(2*i-1:2*i+1,2*i-1:2*i+1)=m(2*i-1:2*i+1,2*i-1:2*i+1)+quadraticReactionElem(1,i,msh);
    % local mass matrix calculated using the reaction solver with lamda = 1    
    k(2*i-1:2*i+1,2*i-1:2*i+1)=k(2*i-1:2*i+1,2*i-1:2*i+1)+QuadraticDiffusionElem(D,i,msh);
    end
    
else error('Invalid element type')
end
    
a0 = m - (1 - theta) * dt * k; % calculates coefficient of c for current timestep
a1 = m + theta * dt * k; % calculates coefficient of c for next timestep

% set boundary conditions for next timestep
a1(1,:) = 0;
a1(end,:) = 0;
a1(1,1) = 1;
a1(end,end) = 1;

for i = 1:tmax/dt % loop accross each timestep
    RHS = a0 * c(:,i); % calculate the right hand side of equation
    RHS(1,1) = c0; % set boundary condition at x = x0
    RHS(end,1) = c1; % set boundary condition at x = x1
    c(:,i+1) = a1\RHS; % calculate solution at next timestep
    
    if i*dt == 0.05 || i*dt == 0.1 || i*dt == 0.3 || i*dt == 1
        plot(x,c(:,i+1))
        hold on
    end
end
hold off
ylabel({'Concentration'});

% Create xlabel
xlabel({'Distance (m)'});

grid on

legend('0.05s','0.1s','0.3s','1s')