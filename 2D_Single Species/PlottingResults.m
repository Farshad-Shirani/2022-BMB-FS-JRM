clear all
close all

load('Results\sol_heterogeneousCapacity.mat')

colorMap1 = [linspace(160/256,256/256,64)' linspace(50/256,210/256,64)' linspace(0/256,145/256,64)' ]; % orange
colorMap2 = [linspace(75/256,220/256,64)' linspace(110/256,256/256,64)' linspace(0/256,160/256,64)' ];  % green


x1_0 = simulationParameters.x1_0;
x1_I = simulationParameters.x1_I;
x2_0 = simulationParameters.x2_0;
x2_J = simulationParameters.x2_J;
Dx1 = discretizationParamaters.Dx1;
Dx2 = discretizationParamaters.Dx2;

x1 = x1_0 : Dx1 : x1_I;
x2 = x2_0 : Dx2 : x2_J;
I = length(x1);
J = length(x2);

[X1, X2] = meshgrid(x1(1:I), x2(1:J-1));



figure, fig_pop1_density = axes;
n1_plotHandle = surf(fig_pop1_density, X1, X2, population.density(:,:,end)');
xlabel(fig_pop1_density, '$x_1$','Interpreter','latex','FontSize', 12);
ylabel(fig_pop1_density, '$x_2$','Interpreter','latex','FontSize', 12);
zlabel(fig_pop1_density, '$n_1$','Interpreter','latex','FontSize', 12);
%set(fig_pop1_density,'ZLim', [0,12]);
colormap(fig_pop1_density, colorMap1);

figure, fig_pop1_trait_mean = axes;
q1_plotHandle = surf(fig_pop1_trait_mean, X1, X2, population.trait_mean(:,:,end)');
xlabel(fig_pop1_trait_mean, '$x_1$','Interpreter','latex','FontSize', 12);
ylabel(fig_pop1_trait_mean, '$x_2$','Interpreter','latex','FontSize', 12);
zlabel(fig_pop1_trait_mean, '$q_1$','Interpreter','latex','FontSize', 12);
%set(fig_pop1_trait_mean,'ZLim', [0,12]);
colormap(fig_pop1_trait_mean, colorMap1);

figure, fig_pop1_trait_variance = axes;
v1_plotHandle = surf(fig_pop1_trait_variance, X1, X2, population.trait_variance(:,:,end)');
xlabel(fig_pop1_trait_mean, '$x_1$','Interpreter','latex','FontSize', 12);
ylabel(fig_pop1_trait_mean, '$x_2$','Interpreter','latex','FontSize', 12);
zlabel(fig_pop1_trait_mean, '$v_1$','Interpreter','latex','FontSize', 12);
set(fig_pop1_trait_variance,'ZLim', [0,5]);
colormap(fig_pop1_trait_variance, colorMap1);
