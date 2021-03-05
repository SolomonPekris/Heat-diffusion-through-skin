function[gamma] = TissueDamage(flow, c0)
%%TISSUEDAMAGE calculates the level of tissue damage that will occure for
%%a given temperature at the outer skin
%
%   [GAMMA] = TISSUEDAMAGE(FLOW,C0) calculates tissue damage GAMMA at the
%   skin boundaries of the dermis, GAMMA.DERM, and at the epidermis,
%   GAMMA.EPI. Blood FLOW can be included using 'yes' or 'no', and the
%   outer temperature is given by C0. 

% set constants
dt = 0.01; % set timestep
Ne6 = 5; % set number of elements divided by 6
c = SkinTempDiffusion1D(dt,flow,Ne6,c0,'no'); % calculates transient temperature distribution
e = 21; % node at epidermis
d = 61; % node at dermis
n = 0; % initialise timestep for epidermis
m = 0; % initialise timestep for dermis
gamma.epi = 0; % initialise gamma
gamma.derm = 0; % initialise derm

for i = 1:length(c) % loop over each temperature value
    if c(e,i) >= 317.15 % check if temp is high enough to cause damage at epidermis
        n = n + 1; 
        gamma.epi(n) = 2 * 10^98 * exp(-12017/(c(e,i)-273.15)); % calculate each point for trapezium
    end
    if c(d,i) >= 317.15 % check if temp is high enough to cause damage at dermis
        m = m + 1; 
        gamma.derm(m) = 2 * 10^98 * exp(-12017/(c(d,i)-273.15)); % calculate each point for trapezium
    end
end

gamma.epi = trapz(gamma.epi) * dt; % numerically integrate to calculate gamma at epidermis
gamma.derm = trapz(gamma.derm) * dt; % numerically integrate to calculate gamma at dermis
    