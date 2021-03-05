%%NODAMAGEBOUNDARYTEMP determines the maximum boundary temperature at x = 0
%%that will not lead to second or third degree burns over 50 seconds.
%
%   NODAMAGEBOUNDARYTEMP checks for what boundary
%   temperature at the skin gamma is 1 or lower at the boundary for the
%   epidermis and the dermis when the temperature is set for 50 seconds. At
%   this value or lower in the epidermis, second-degree burns will not
%   occur and at the dermis, third-degree burns will not occur.
clear
m = 0; % initialises stepping variable
T = 317.15:2:393.15; % sets temperature values

for c0 = T % loops through all potentially burning temperatures
    G = TissueDamage("no", c0); % does not include flow
    m = m + 1;
    gamma.epi(m) = G.epi; % inputs current iterations solution for epidermis
    gamma.derm(m) = G.derm; % inputs current iterations solution for dermis
end

m = 0;

for c0 = T
    G = TissueDamage("yes", c0); % same as previous loop but includes bloodflow
    m = m + 1;
    gammaflow.epi(m) = G.epi;
    gammaflow.derm(m) = G.derm;
end

plot(T,gamma.derm,T,gammaflow.derm,'color',[0.00,0.45,0.74],'LineWidth',1) % plots dermis results
hold on
plot(T,gamma.epi,T,gammaflow.epi,'color',[1 0 0],'LineWidth',1) % plots epidermis results
plot([393.15 310],[1 1],'color',[0 0 0], 'LineWidth',1) % plots line of gamma=1

% Create ylabel
ylabel({'Gamma'});

% Create xlabel
xlabel({'Boundary Temperature (K)'});

%LineWidth('1');

grid on

legend('Dermis no flow','Dermis flow','Epidermis no flow','Epidermis flow','Gamma = 1')

set(gca, 'XScale', 'log', 'YScale', 'log'); % set to log scale
hold off

