function[c] = SkinTempDiffusion1D(dt,flow,Ne6,c0,plt)
%%SKINTEMPDIFFUSION1D calculated the transient temperature response of
%%different layers of skin in 1D using FEM.
%
%   [C] = SKINTEMPDIFFUSION1D(DT,FLOW,NE6,C0,PLT) calculates the temperature C
%   at different points in the skin over time when given a step input at
%   the outer boundary C0 by modelling it as a transient diffusion reaction 
%   equation with a source term. It is then solved as an FEM.
%   The backward Euler stepping method is used with
%   timestep DT. FLOW determines if blood flow is considered using the
%   inputs 'yes' or 'no'. The number of elements must be divisible by 6 so
%   that a node exists at each of the skin boundary layers, so NE6 is the
%   number of elements to be used divided by 6. Quadratic basis functions
%   are used to calculate the stiffness and mass matrices, and the source
%   vector. PLT determins if a plot will be made.

% set constants. suffix epi is epidermis, derm is dermis, and subq is
% sub-cutaneous
Ne = Ne6 * 6; % multiply number of elements by six so that there is a node at each skin boundary layer
ct0 = 310.15; % initial temperature of mesh
c1 = 310.15; % dirichlet boundary at x = xb

% thermal conductivities
kepi = 25;
kderm = 40;
ksubq = 20;

% blood flow
if flow == 'yes'
    g = 0.0375;
elseif flow == 'no'
    g = 0;
else 
    error('Invalid flow condition')
end

ro = 1200; % tissue density for all layers
cp = 3300; % specific heat capacity for all tissue layers
rob = 1060; % blood density
cpb = 3770; % specific heat capacity for blood
tb = 310.15; % blood temperature

tmax = 50; % set end time
c = zeros(2*Ne+1,(1+tmax/dt)); % initialise solution matrix
x = (0:0.01/(2*Ne):0.01)'; % initialise distance vector
c(:,1) = 310.15; % set initial temperature conditions
theta = 0.5; % select time-stepping method. Theta = 0.5 is Crank-Nicolson
m = zeros(2*Ne+1,2*Ne+1); % Initialise global mass matrix
k = zeros(2*Ne+1,2*Ne+1); % Initiate global diffusion matrix
s = zeros(2*Ne+1,1); % Initiate global source vector
msh = OneDimLinearMeshGen(0,0.01,Ne);

for i = 1:Ne % Loop through each element calculating the global mass matrix, global stiffness matrix and global source vector
    if i <= Ne/6 % sets conditions at epidermis
        D = kepi / (ro * cp); % set diffusion coefficient
        lamda = 0; % set lamda
        f = 0; % set source term
    elseif (i <= Ne/2) && (i > Ne/6) % set tissue conditions at dermis
        D = kderm / (ro * cp); % set diffusion coefficient
        f = g * rob * cpb * tb/ (ro * cp); % set source term
        lamda = -f / tb; % set lamda
    else % mesh is at subcutaneous set material parameters
        D = ksubq / (ro * cp); % set diffusion coefficient
        f = g * rob * cpb * tb / (ro * cp); % set source term
        lamda = -f / tb; % set lamda
    end
    % local mass matrix calculated using the reaction solver with lamda = 1
    m(2*i-1:2*i+1,2*i-1:2*i+1) = m(2*i-1:2*i+1,2*i-1:2*i+1) + quadraticReactionElem(1,i,msh);
    % adds local diffusion and reaction matrices to global stifness matrix
    k(2*i-1:2*i+1,2*i-1:2*i+1) = k(2*i-1:2*i+1,2*i-1:2*i+1) + QuadraticDiffusionElem(D,i,msh) - quadraticReactionElem(lamda,i,msh);
    % adds local source term to global source vector
    s(2*i-1:2*i+1) = s(2*i-1:2*i+1) + QuadraticSource(f,i,msh);
end

a0 = m - (1 - theta) * dt * k; % calculates coefficient of c for current timestep
a1 = m + theta * dt * k; % calculates coefficient of c for next timestep

% set boundary conditions for next timestep
a1(1,:) = 0;
a1(end,:) = 0;
a1(1,1) = 1;
a1(end,end) = 1;

for i=1:tmax/dt % loop over time steps
    RHS = a0 * c(:,i) + dt * s; % calculate the right hand side of equation
    RHS(1,1) = c0; % set boundary condition at x = x0
    RHS(end,1) = c1; % set boundary condition at x = x1
    c(:,i+1) = a1\RHS; % calculate solution at next timestep
    if plt == "yes"
        if mod(log2(i),1) == 0 || i == tmax/dt
            plot(x,c(:,i))
            hold on
        end
    elseif plt == "no"
        continue
    else
        error('Invalid plot input')
    end
end

if plt == "yes"
    xlabel('Distance (m)')
    ylabel('Temperature (K)')

    hold off
end
