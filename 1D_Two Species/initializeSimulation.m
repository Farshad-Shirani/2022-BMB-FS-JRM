function [initialPopulations, modelParameters, simulationParameters, discretizationParamaters, epsilon, eta] = initializeSimulation

epsilon = 0;%1e-23; % threshol for smallest value of population density (to avoid singularity due to 1/N term) 

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

%---model parameters---------------------------------------------------
global D1 D2 V_s V_u1 V_u2 U kappa K1 K2 theta Q_opt R1 R2
D1 = 1;  % diffusion matrix for populationn 1
D2 = 1;  % diffusion matrix for populationn 2
V_s = 1/0.2;    % 1/measure of the strength of stabilizing selection (V_s = 1/S, where S defined in Table 1)
V_u1 = 4;    % variance of the within phenotipic-resource utility curve for population 1 (V_u1 = V_1 where V is defined in Table 1)
V_u2 = 4;    % variance of the within phenotipic-resource utility curve for population 2 (V_u2 = V_2 where V is defined in Table 1)
U = 0;  % rate of increase in trait variance due to mutation
kappa = 0;  % population impact factor
K1 = 1 * ones(I, 1);  % carrying capacity for population 1
K2 = 1 * ones(I, 1);  % carrying capacity for population 2
R1 = 2 * ones(I, 1);  % intrinsic rate of increase for optimally adapted individuals of population 1
R2 = 2 * ones(I, 1);  % intrinsic rate of increase for optimally adapted individuals of population 2
theta = 0.2; % optimum trait gradient
Q_0 = 10; % optimum trait at x_0
Q_opt = Q_0 + theta  * (x - x(1))';    % optimum trait value  
%Q_opt = Q_0 + theta  * (x - x(1))' + ( (10-theta)/2 * log(1+exp(2*(x - 45)  )) )';    % optimum trait value (sharp change to dQ_opt = 10 near boundary )
modelParameters = struct('D1', D1, 'D2', D2,  'V_s', V_s, 'V_u1', V_u1, 'V_u2', V_u2, 'U', U, 'kappa', kappa, 'K1', K1, 'K2', K2, 'Q_opt', Q_opt, 'R1', R1, 'R2', R2);

%---solver parameters----------------------------------------------------
eta = 1/2;

%---initial values-------------------------------------------------------
initialPopulations = struct('density', [], 'trait_mean', [], 'trait_variance', []);

%---For calculating nominal solution------------------------------------------------------------ 
center = -15; % population center
radius = sqrt(2 * 1); % population radius in each direction 
initialPopulations(1).density = 0.5 * sech( abs( x - center ) / radius )' + epsilon;
center_index_x = dsearchn(x', center); 
initialPopulations(1).trait_mean = 0.6 * ( Q_opt - Q_opt(center_index_x) ) + Q_opt(center_index_x);
initialPopulations(1).trait_variance = 1 * ones(I, 1);

center = 15; % population center 
radius = sqrt(2 * 1); % population radius in each direction
initialPopulations(2).density = 0.5 * sech( abs( x - center ) / radius )' + epsilon;
center_index_x = dsearchn(x', center); 
initialPopulations(2).trait_mean = 0.6 * ( Q_opt - Q_opt(center_index_x) ) + Q_opt(center_index_x);
initialPopulations(2).trait_variance = 1 * ones(I, 1);

% %---For introduction of a new species to already established species --------------------------- 
% load('Results\sol_single_nominal.mat', 'population' ) 
% initialPopulations(1).density = population.density(:,end);
% initialPopulations(1).trait_mean = population.trait_mean(:,end);
% initialPopulations(1).trait_variance = population.trait_variance(:,end);
% clear population
% 
% center = 0; % population center
% radius = sqrt(2 * 1); % population radius in each direction
% initialPopulations(2).density = 0.5 * sech( abs( x - center ) / radius )' + epsilon;
% center_index_x = dsearchn(x', center); 
% initialPopulations(2).trait_mean = 0.6 * ( Q_opt - Q_opt(center_index_x) ) + Q_opt(center_index_x);
% initialPopulations(2).trait_variance = 1 * ones(I, 1);

% %---For stable range/extinction near sharp Q_opt-------------------------------------------- 
% center = 10; % population center
% radius = sqrt(2); % population radius in each direction 
% initialPopulations(1).density = 0.5 * sech( abs( x - center ) / radius )' + epsilon;
% center_index_x = dsearchn(x', center); 
% initialPopulations(1).trait_mean = 1 * ( Q_opt - Q_opt(center_index_x) ) + Q_opt(center_index_x);
% initialPopulations(1).trait_variance = 1 * ones(I, 1);
% 
% center = 35; % population center 
% radius = sqrt(2); % population radius in each direction
% initialPopulations(2).density = 0.5 * sech( abs( x - center ) / radius )' + epsilon;
% center_index_x = dsearchn(x', center); 
% initialPopulations(2).trait_mean = 1 * ( Q_opt - Q_opt(center_index_x) ) + Q_opt(center_index_x);
% initialPopulations(2).trait_variance = 1 * ones(I, 1);


end

