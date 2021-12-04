function [initialPopulation, modelParameters, simulationParameters, discretizationParamaters, epsilon, eta] = initializeSimulation

epsilon = 0;%1e-16; % threshol for smallest value of population density (to avoid singularity due to 1/N term) 

%---simulation parameters----------------------------------------------
T = 50; %final time
storageTimes = 0 : 0.1 : T; % solution storage times. Solutins are sored at closest time samples depending on Dt
x_0 = -50;    x_I = 50; % x1 range
simulationParameters = struct('T', T, 'times', storageTimes, 'x_0', x_0, 'x_I', x_I);

%---discretization parameters------------------------------------------
Dt = 0.002;  % time steps
Dx = (x_I - x_0)/ 1000; %x1 mesh   
discretizationParamaters = struct('Dt', Dt, 'Dx', Dx);

x = x_0 : Dx : x_I;
I = length(x);

%---model parameters----------------------------------------------------
global D V_s V_u U kappa K theta Q_opt R 
D = 1;  % diffusion matrix 
V_s = 1/0.2;    % 1/measure of the strength of stabilizing selection (V_s = 1/S, where S defined in Table 1)
V_u = 4;    % variance of the within phenotipic-resource utility curve (V_u = V where V is defined in Table 1)
U = 0;  % rate of increase in trait variance due to mutation
kappa = 0;  % population impact factor
K = 1 * ones(I, 1);  % constant carrying capacity
%K = 0.01 - (x - x_0).*(x - x_I) ./ ( (x_I - x_0)/2 ).^2;  % parabolic carrying capacity
%K = 1.1 + 0.01 +(x - x_0).*(x - x_I) ./ ( (x_I - x_0)/2 ).^2;  % parabolic carrying capacity
%K = 0.55 + 0.45 * sin( 16*pi/(x_I) * x );
R = 2 * ones(I, 1);  % intrinsic rate of increase for optimally adapted individuals
theta = 0.2; % optimum trait gradient dQ
Q_0 = 10; % optimum trait at x_0
Q_opt = Q_0 + theta  * (x - x(1))';    % optimum trait value  
modelParameters = struct('D', D, 'V_s', V_s, 'V_u', V_u, 'U', U, 'kappa', kappa, 'K', K, 'Q_opt', Q_opt, 'R', R);

%---solver parameters----------------------------------------------------
eta = 1/2;

%% Initial values --------------------------------------------------------
initialPopulation = struct('density', [], 'trait_mean', [], 'trait_variance', []);
center = 0; % population center

%---For calculating nominal solution------------------------------------------------------------ 
radius = sqrt(2); % population radius in each direction 
initialPopulation.density = 0.5 * sech( abs(x - center) / radius )' + epsilon;
center_index_x = dsearchn(x', center); 
initialPopulation.trait_mean = 0.6 * ( Q_opt - Q_opt(center_index_x) ) + Q_opt(center_index_x);
initialPopulation.trait_variance = 1 * ones(I, 1);

%---For calculationg solution at different environmental gradients-------------
% radius = sqrt(2); % population radius in each direction 
% initialPopulation.density = 1 * sech( abs(x - center) / radius )' + epsilon;
% center_index_x = dsearchn(x', center); 
% initialPopulation.trait_mean = 1 * ( Q_opt - Q_opt(center_index_x) ) + Q_opt(center_index_x);
% initialPopulation.trait_variance = 1 * ones(I, 1);
 

%---Initialize with the nominal solution at t=4
% load('Results\sol_nominal.mat', 'population' ) 
% %load('Results\sol_gradient_4.mat', 'population' ) 
% initialPopulation.density = population.density(:,41);
% initialPopulation.trait_mean = population.trait_mean(:,41);
% initialPopulation.trait_variance = population.trait_variance(:,41);
% clear population

 
% %---Initialize with the parabollic carrying capacity solution at t=50
% load('Results\sol_ParabolicCapacity.mat', 'population' ) 
% initialPopulation.density = population.density(:,end);
% initialPopulation.trait_mean = population.trait_mean(:,end);
% initialPopulation.trait_variance = population.trait_variance(:,end);
% clear population




