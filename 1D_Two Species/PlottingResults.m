clear all
close all

load('Results\sol_nominal.mat')

% Colors
darkBlue =  [0, 114/255, 189/255];
lightBlue = [0.3010, 0.7450, 0.9330];
transparentBlue =  [189/255, 223/255, 246/255];
orange = [0.8500, 0.3250, 0.0980];
yellow =  [0.9290, 0.6940, 0.1250];
purple =  [0.4940, 0.1840, 0.5560];
green =  [0.4660, 0.6740, 0.1880];
darkGreen =  [33/255, 186/255, 140/255];
transparentGreen = [191/255, 242/255, 226/255];
darkRed = [162/255, 20/255, 47/255];
transparentRed = [247/255, 210/255, 217/255];
colors = [darkBlue; orange; yellow; purple; green; lightBlue; darkRed ];

colorMap1 = [linspace(160/256,256/256,64)' linspace(50/256,210/256,64)' linspace(0/256,145/256,64)' ]; % orange
colorMap2 = [linspace(75/256,220/256,64)' linspace(110/256,256/256,64)' linspace(0/256,160/256,64)' ];  % green


x_0 = simulationParameters.x_0;
x_I = simulationParameters.x_I;
Dx = discretizationParamaters.Dx;

x = x_0 : Dx : x_I;
I = length(x);
numSamples = length(simulationParameters.times);

N1 = populations(1).density;
Q1 = populations(1).trait_mean;
V1 = populations(1).trait_variance;
N2 = populations(2).density;
Q2 = populations(2).trait_mean;
V2 = populations(2).trait_variance;

figure, 
plot(x, N1(:,[1:floor(numSamples/50):numSamples])','Color', orange, 'LineWidth', 0.5);
hold on
plot(x, N2(:,[1:floor(numSamples/50):numSamples])','Color', darkGreen, 'LineWidth', 0.5);
plot(x, N2(:,1)','Color', darkGreen, 'LineWidth', 1.5);
plot(x, N1(:,1)','Color', orange, 'LineWidth', 1.5);
plot(x, N2(:,81)','Color', darkBlue, 'LineWidth', 2);
plot(x, N1(:,81)','Color', darkRed, 'LineWidth', 2);
hold off
xlabel('Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel('Population Density $[\mathtt{N}/\mathtt{X}]$','Interpreter','latex','FontSize', 12);

figure, 
plot(x, Q1(:,[1:floor(numSamples/25):numSamples])','Color', orange, 'LineWidth', 0.5);
hold on
plot(x, Q2(:,[1:floor(numSamples/25):numSamples])','Color', darkGreen, 'LineWidth', 0.5);
plot(x, Q2(:,1)','Color', darkGreen, 'LineWidth', 1.5);
plot(x, Q1(:,1)','Color', orange, 'LineWidth', 1.5);
plot(x, Q2(:,81)', 'Color', darkBlue, 'LineWidth', 2);
plot(x, Q1(:,81)', 'Color', darkRed, 'LineWidth', 2);
plot(x, modelParameters.Q_opt','k', 'LineWidth', 1);
hold off
xlabel('Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel('Trait Mean $[\mathtt{Q}]$','Interpreter','latex','FontSize', 12);

figure, 
plot(x, V1(:,[1:floor(numSamples/25):numSamples])', 'Color', orange, 'LineWidth', 0.5);
hold on
plot(x, V2(:,[1:floor(numSamples/25):numSamples])', 'Color', darkGreen, 'LineWidth', 0.5);
% plot(x, V2(:,1)','Color', darkGreen, 'LineWidth', 1.5);
plot(x, V1(:,1)','Color', orange, 'LineWidth', 1.5);
plot(x, V1(:,81)','Color', darkRed, 'LineWidth', 2);
hold off
xlabel('Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel('Trait Variance $[\mathtt{Q}^2]$','Interpreter','latex','FontSize', 12);

figure, 
plot(x, Q1(:,end)','Color', transparentRed, 'LineWidth', 1.5);
hold on
plot(x, Q2(:,end)','Color', transparentBlue, 'LineWidth', 1.5);
Q1_edgeIndex = dsearchn(x', 19); % N1(:,end)=0 at approximately x=20
Q2_edgeIndex = dsearchn(x', -19); % N2(:,end)=0 at approximately x=-20
plot(x(1:Q1_edgeIndex), Q1(1:Q1_edgeIndex,end)','Color', darkRed, 'LineWidth', 1.5);
plot(x(Q2_edgeIndex:end), Q2(Q2_edgeIndex:end,end)','Color', darkBlue, 'LineWidth', 1.5);
plot(x, modelParameters.Q_opt','k', 'LineWidth', 1);
hold off
grid on
xlabel('Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel('Trait Mean $[\mathtt{Q}]$','Interpreter','latex','FontSize', 12);



