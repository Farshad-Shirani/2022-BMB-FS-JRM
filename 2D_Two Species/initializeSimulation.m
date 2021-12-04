function [initialPopulations, modelParameters, simulationParameters, discretizationParamaters, epsilon, eta] = initializeSimulation

epsilon = 1e-14; % threshol for smallest value of population density (to avoid singularity due to 1/N term) 
%epsilon = 0; % Works for the simulation presented in the paper (introduction of a new species to already established species)

%---simulation parameters----------------------------------------------
T = 50; %final time
storageTimes = 0 : 5 : T; % solution storage times. Solutins are sored at closest time samples depending on Dt
x1_0 = -50;    x1_I = 50; % x1 range
x2_0 = -30;    x2_J = 30; % x2 range
%x2_0 = -50;    x2_J = 50; % x2 range
simulationParameters = struct('T', T, 'times', storageTimes, 'x1_0', x1_0, 'x1_I', x1_I, 'x2_0', x2_0, 'x2_J', x2_J);

%---discretization parameters------------------------------------------
Dt = 0.01;  % time steps
Dx1 = (x1_I - x1_0)/ 200; %x1 mesh   
Dx2 = (x2_J - x2_0)/ 120; %x2 mesh 
%Dx2 = (x2_J - x2_0)/ 200; %x2 mesh  
discretizationParamaters = struct('Dt', Dt, 'Dx1', Dx1, 'Dx2', Dx2);

x1 = x1_0 : Dx1 : x1_I;
x2 = x2_0 : Dx2 : x2_J;
I = length(x1);
J = length(x2);

%---model parameters----------------------------------------------------
global D1 D2 V_s V_u1 V_u2 U kappa K1 K2 theta Q_opt R1 R2
D1 = diag([1 1]);  % diffusion matrix for populationn 1
D2 = diag([1 1]);  % diffusion matrix for populationn 2
% D2 = 6 * diag([1 1]);  % diffusion matrix for populationn 2
V_s = 1/0.2;    % 1/measure of the strength of stabilizing selection (V_s = 1/S, where S defined in Table 1)
V_u1 = 4;    % variance of the within phenotipic-resource utility curve for population 1 (V_u1 = V_1 where V is defined in Table 1)
V_u2 = 4;     % variance of the within phenotipic-resource utility curve for population 2 (V_u2 = V_2 where V is defined in Table 1)
%V_u2 = 6;     % variance of the within phenotipic-resource utility curve for population 2 (V_u2 = V_2 where V is defined in Table 1)
U = 0;  % rate of increase in trait variance due to mutation
kappa = 0;  % population impact factor
K1 = 1 * ones(I, J-1);  % carrying capacity for population 1
% [X1, X2] = meshgrid(x1(1:I), x2(1:J-1));
% K1 = ( 1 - 0.5 * (1 - exp(-(X1+50)/30 ) )  .* exp( 1 ./ ( abs(X2/50) - 1 -eps ) ) .*  (1+ sin(pi*sqrt( X2.^2/1000 + (X1+80).*(X2+50) ) / 4)) )';
K2 = 1 * ones(I, J-1);  % carrying capacity for population 2
R1 = 2 * ones(I, J-1);  % intrinsic rate of increase for optimally adapted individuals of population 1
R2 = 2 * ones(I, J-1);  % intrinsic rate of increase for optimally adapted individuals of population 2
theta = 0.2; % optimum trait gradient
Q_0 = 10; % optimum trait at x1_0
Q_opt = Q_0 + theta  * (x1 - x1(1))' * ones(1, J-1);    % optimum trait value  
modelParameters = struct('D1', D1, 'D2', D2,  'V_s', V_s, 'V_u1', V_u1, 'V_u2', V_u2, 'U', U , 'kappa', kappa, 'K1', K1, 'K2', K2, 'Q_opt', Q_opt, 'R1', R1, 'R2', R2);

%---solver parameters----------------------------------------------------
eta = 1/2;

%---initial values--------------------------------------------------------
initialPopulations = struct('density', [], 'trait_mean', [], 'trait_variance', []);
[X1, X2] = meshgrid(x1(1:I), x2(1:J-1));

%---For calculating nominal solution------------------------------------------------------------ 
center = [-15 0]; % population center
radius = 2 * [1; 1]; % population radius in each direction 
initialPopulations(1).density = 0.5 * sech( sqrt( ( X1 - center(1) ).^2 / radius(1) + ( X2 - center(2) ).^2 / radius(2) ) )' + epsilon;
center_index_x1 = dsearchn(x1', center(1)); 
center_index_x2 = dsearchn(x2', center(2));
initialPopulations(1).trait_mean = 0.6 * ( Q_opt - Q_opt(center_index_x1, center_index_x2) ) + Q_opt(center_index_x1, center_index_x2);
initialPopulations(1).trait_variance = 1 * ones(I, J-1);

center = [15 0]; % population center 
radius = 2 * [1; 1]; % population radius in each direction
initialPopulations(2).density = 0.5 * sech( sqrt( ( X1 - center(1) ).^2 / radius(1) + ( X2 - center(2) ).^2 / radius(2) ) )' + epsilon;
center_index_x1 = dsearchn(x1', center(1)); 
center_index_x2 = dsearchn(x2', center(2));
initialPopulations(2).trait_mean = 0.6 * ( Q_opt - Q_opt(center_index_x1, center_index_x2) ) + Q_opt(center_index_x1, center_index_x2);
initialPopulations(2).trait_variance = 1 * ones(I, J-1);

% %---For introduction of a new species to already established species --------------------------- 
% load('Results\sol_single_heterogeneousCapacity.mat', 'population' ) 
% initialPopulations(1).density = population.density(:,:,end);
% initialPopulations(1).trait_mean = population.trait_mean(:,:,end);
% initialPopulations(1).trait_variance = population.trait_variance(:,:,end);
% clear population
% 
% center = [0 0]; % population center 
% radius = 2 * [1; 1]; % population radius in each direction
% %radius = [20 200]; % population radius in each direction
% initialPopulations(2).density = 0.1 * sech( sqrt( ( X1 - center(1) ).^2 / radius(1) + ( X2 - center(2) ).^2 / radius(2) ) )' + epsilon;
% center_index_x1 = dsearchn(x1', center(1)); 
% center_index_x2 = dsearchn(x2', center(2));
% initialPopulations(2).trait_mean = 0.6 * ( Q_opt - Q_opt(center_index_x1, center_index_x2) ) + Q_opt(center_index_x1, center_index_x2);
% initialPopulations(2).trait_variance = 1 * ones(I, J-1);

end

