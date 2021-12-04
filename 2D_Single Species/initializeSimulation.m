function [initialPopulation, modelParameters, simulationParameters, discretizationParamaters, epsilon, eta] = initializeSimulation

epsilon = 1e-16; % threshol for smallest value of population density (to avoid singularity due to 1/N term) 

%---simulation parameters----------------------------------------------
T = 100; %final time
storageTimes = 0 : 5 : T; % solution storage times. Solutins are sored at closest time samples depending on Dt
x1_0 = -50;    x1_I = 50; % x1 range
x2_0 = -50;    x2_J = 50; % x2 range
simulationParameters = struct('T', T, 'times', storageTimes, 'x1_0', x1_0, 'x1_I', x1_I, 'x2_0', x2_0, 'x2_J', x2_J);

%---discretization parameters------------------------------------------
Dt = 0.01;  % time steps
Dx1 = (x1_I - x1_0)/ 200; %x1 mesh   
Dx2 = (x2_J - x2_0)/ 200; %x2 mesh  
discretizationParamaters = struct('Dt', Dt, 'Dx1', Dx1, 'Dx2', Dx2);

x1 = x1_0 : Dx1 : x1_I;
x2 = x2_0 : Dx2 : x2_J;
I = length(x1);
J = length(x2);

%---model parameters----------------------------------------------------
global D V_s V_u U kappa K theta Q_opt R 
D = diag([1 1]);  % diffusion matrix 
V_s = 1/0.2;    % 1/measure of the strength of stabilizing selection (V_s = 1/S, where S defined in Table 1)
V_u = 4;    % variance of the within phenotipic-resource utility curve (V_u = V where V is defined in Table 1)
U = 0;  % rate of increase in trait variance due to mutation
kappa = 0;  % population impact factor
% K = 1 * ones(I, J-1);  % carrying capacity
[X1, X2] = meshgrid(x1(1:I), x2(1:J-1));
K = ( 1 - 0.5 * (1 - exp(-(X1+50)/30 ) )  .* exp( 1 ./ ( abs(X2/50) - 1 -eps ) ) .*  (1+ sin(pi*sqrt( X2.^2/1000 + (X1+80).*(X2+50) ) / 4)) )';
R = 2 * ones(I, J-1);  % intrinsic rate of increase for optimally adapted individuals
theta = 0.2; % optimum trait gradient
Q_0 = 10; % optimum trait at x1_0
Q_opt = Q_0 + theta  * (x1 - x1(1))' * ones(1, J-1);    % optimum trait value  
modelParameters = struct('D', D, 'V_s', V_s, 'V_u', V_u, 'U', U, 'kappa', kappa, 'K', K, 'Q_opt', Q_opt, 'R', R);

%---solver parameters----------------------------------------------------
eta = 1/2;

%---initial values--------------------------------------------------------
initialPopulation = struct('density', [], 'trait_mean', [], 'trait_variance', []);
[X1, X2] = meshgrid(x1(1:I), x2(1:J-1));

center = [0 0]; % population center
radius = 2 * [1; 1]; % population radius in each direction 
initialPopulation.density = 1 * sech( sqrt( ( X1 - center(1) ).^2 / radius(1) + ( X2 - center(2) ).^2 / radius(2) ) )' + epsilon;

center_index_x1 = dsearchn(x1', center(1)); 
center_index_x2 = dsearchn(x2', center(2));
initialPopulation.trait_mean = 0.9 * ( Q_opt - Q_opt(center_index_x1, center_index_x2) ) + Q_opt(center_index_x1, center_index_x2);

initialPopulation.trait_variance = 1 * ones(I, J-1);

end

