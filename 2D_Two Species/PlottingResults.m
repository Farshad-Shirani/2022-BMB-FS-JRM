clear all
close all

load('Results\sol_heterogeneous_capacity.mat')

% colorMap1 = [linspace(188/255,255/255,64)' linspace(56/255,231/255,64)' linspace(5/256,209/255,64)' ]; % orange
% colorMap2 = [linspace(8/255,196/255,64)' linspace(165/256,252/255,64)' linspace(118/256,235/255,64)' ];  % green
%---Reverse color maps
colorMap1 = [linspace(255/255, 188/255, 64)' linspace(231/255, 56/255, 64)' linspace(209/255, 5/256, 64)' ]; % orange
colorMap2 = [linspace(196/255, 8/255, 64)' linspace(252/255, 165/256, 64)' linspace(235/255, 118/256, 64)' ];  % green

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

N1 = populations(1).density;
Q1 = populations(1).trait_mean;
V1 = populations(1).trait_variance;
N2 = populations(2).density;
Q2 = populations(2).trait_mean;
V2 = populations(2).trait_variance;


%% Plot sample frames ==================================================================
k = 41; % frame number

figure, fig_pop1_density = axes;
n1_plotHandle = surface(fig_pop1_density, X1, X2, N1(:,:,k)');
% xlabel(fig_pop1_density, '$x_1$','Interpreter','latex','FontSize', 12);
% ylabel(fig_pop1_density, '$x_2$','Interpreter','latex','FontSize', 12);
set(n1_plotHandle,'EdgeColor', 'none');
caxis(fig_pop1_density, [min(N1(:)),max(N1(:))] )
colorbar
axis off
axis image
colormap(fig_pop1_density, colorMap1);

figure, fig_pop2_density = axes;
n2_plotHandle = surface(fig_pop2_density, X1, X2, N2(:,:,k)');
% xlabel(fig_pop2_density, '$x_1$','Interpreter','latex','FontSize', 12);
% ylabel(fig_pop2_density, '$x_2$','Interpreter','latex','FontSize', 12);
set(n2_plotHandle,'EdgeColor', 'none');
caxis(fig_pop2_density, [min(N2(:)),max(N2(:))] )
colorbar
axis off
axis image
colormap(fig_pop2_density, colorMap2);

