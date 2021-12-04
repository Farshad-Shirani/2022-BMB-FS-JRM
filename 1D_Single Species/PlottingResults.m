clear all
close all

load('Results\sol_nominal.mat')

% Colors
darkBlue =  [0, 0.4470, 0.7410];
orange = [0.8500, 0.3250, 0.0980];
yellow =  [0.9290, 0.6940, 0.1250];
purple =  [0.4940, 0.1840, 0.5560];
green =  [0.4660, 0.6740, 0.1880];
darkGreen =  [33/255, 186/255, 140/255];
lightBlue = [0.3010, 0.7450, 0.9330];
darkRed = [0.6350, 0.0780, 0.1840];
colors = [darkBlue; orange; yellow; purple; green; lightBlue; darkRed ];

colorMap1 = [linspace(160/256,256/256,64)' linspace(50/256,210/256,64)' linspace(0/256,145/256,64)' ]; % orange
colorMap2 = [linspace(75/256,220/256,64)' linspace(110/256,256/256,64)' linspace(0/256,160/256,64)' ];  % green

x_0 = simulationParameters.x_0;
x_I = simulationParameters.x_I;
Dx = discretizationParamaters.Dx;

x = x_0 : Dx : x_I;
I = length(x);
numSamples = length(simulationParameters.times);

N = population.density;
Q = population.trait_mean;
V = population.trait_variance;

figure, 
plot(x, N(:,[1:floor(numSamples/25):numSamples])','Color', orange, 'LineWidth', 0.5);
hold on
plot(x, N(:,1)','Color', orange, 'LineWidth', 1.5);
plot(x, N(:,81)','Color', darkRed, 'LineWidth', 1.5);
hold off
xlabel('Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel('Population Density $[\mathtt{N}/\mathtt{X}]$','Interpreter','latex','FontSize', 12);


figure, 
plot(x, Q(:,[1:floor(numSamples/25):numSamples])','Color', orange, 'LineWidth', 0.5);
hold on
plot(x, Q(:,1)','Color', orange, 'LineWidth', 1.5);
plot(x, Q(:,81)','Color', darkRed, 'LineWidth', 1.5);
plot(x, modelParameters.Q_opt','k', 'LineWidth', 1);
hold off
xlabel('Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel('Trait Mean $[\mathtt{Q}]$','Interpreter','latex','FontSize', 12);

figure, 
plot(x, V(:,[1:floor(numSamples/25):numSamples])', 'Color', orange, 'LineWidth', 0.5);
hold on
plot(x, V(:,1)','Color', orange, 'LineWidth', 1.5);
plot(x, V(:,81)','Color', darkRed, 'LineWidth', 1.5);
hold off
xlabel('Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel('Trait Variance $[\mathtt{Q}^2]$','Interpreter','latex','FontSize', 12);


%% Plot G (Mean growth rate) ====================================================
V_s = modelParameters.V_s;
V_u = modelParameters.V_u;
kappa = modelParameters.kappa;
K = modelParameters.K;
Q_opt = modelParameters.Q_opt;
R = modelParameters.R;

Qm2kVu = Q - 2 * kappa * V_u;

V_squared = V.^2;
Vp2Vu = V + 2 * V_u;
Vp2Vu_squared = Vp2Vu.^2;
VpVp2Vu = V + Vp2Vu;
VpVp2Vu_squared = VpVp2Vu.^2;

Q_opt_squared = Q_opt(:).^2;
Q_squared = Q.^2;
Q_cubed = Q_squared .* Q;
Q_optQ = Q_opt(:) .* Q;
twoQ_optV = 2 * Q_opt(:) .* V;
Q_optQ_squared = Q_opt(:) .* Q_squared;
Q_opt_squaredQ = Q_opt_squared .* Q;
QmQ_opt = Q - Q_opt(:);
VmQ_squared = V - Q_squared;

C = sqrt(2*V_u) * exp( kappa^2 * V_u ) ./ sqrt(VpVp2Vu);
M = exp( -(2 * kappa * V_u)^2 ./ (2 * VpVp2Vu) );
L = ( V .* Qm2kVu + Vp2Vu .* Q ) ./ VpVp2Vu;
S = V .* Vp2Vu ./ VpVp2Vu + L .* (L - 2*Q);
E = ( twoQ_optV + 2 * Q_optQ_squared - Q_opt_squaredQ - 3 * V .* Q - Q_cubed ) / (2*V_s);
Y = ( twoQ_optV .* Q - 2 * Q_opt(:) .* Q_cubed - Q_opt_squared .* VmQ_squared - 3 * V_squared + Q_cubed .* Q ) / (2*V_s);

R_K =  R(:) ./ K(:);
RN_K = R_K .* N;
RCN_K = RN_K .* C;
RMN_K = RN_K .* M;
RMCN_K = RCN_K .* M;
RMC_K = R_K .* M .* C;

G = R(:) - RMCN_K - ( QmQ_opt.^2 + V ) / (2 * V_s);

figure, 
plot(x, G(:,[1:floor(numSamples/25):numSamples])','Color', orange, 'LineWidth', 0.5);
hold on
plot(x, G(:,1)','Color', orange, 'LineWidth', 1.5);
plot(x, G(:,81)','Color', darkRed, 'LineWidth', 1.5);
hold off
xlabel('Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel('Mean Intrinsic Growth $[1/\mathtt{T}]$','Interpreter','latex','FontSize', 12);


