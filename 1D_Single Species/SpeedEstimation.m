clear all
close all

numCurves = 28;
purple =  [0.4940, 0.1840, 0.5560];

invasionSpeed = zeros(numCurves,1);
gradient = zeros(numCurves,1);
amplitude = zeros(numCurves,1);
maxVariance = zeros(numCurves,1);
for i = 1:numCurves
    path = strcat('Results\sol_gradient_', num2str(i) ,'.mat');
    load(path);
    n = population.density;
    amplitude(i) = max(n(:,end));
    gradient(i) = (modelParameters.Q_opt(end) - modelParameters.Q_opt(1)) / (simulationParameters.x_I-simulationParameters.x_0);
    times = simulationParameters.times;
    x = simulationParameters.x_0 : discretizationParamaters.Dx : simulationParameters.x_I;

    n = n(floor(size(n,1)/2):end, :); % only keep half of the symmetric density
    x = x(end-size(n,1)+1:end);
    
    threshold = 0.25 * amplitude(i);
    edgePoints = zeros(1,length(times));
    
    for j = 1:length(times)
        edge_index_x = dsearchn(n(:,j), threshold);
        edgePoints(j)= x(edge_index_x);
    end
    
    noThresholdCrossing = (edgePoints == simulationParameters.x_0);
    edgePoints(noThresholdCrossing) = []; % remove no threshold crossings
    times(noThresholdCrossing) = [];
    
    len = length(edgePoints);
    edgePoints(1:floor(len/5)) = []; %remove transient 20 percent at the begining
    times(1:floor(len/5)) = [];
    edgePoints(end-floor(len/10) : end) = []; % remove nonlinear 10 percent at the end
    times(end-floor(len/10) : end) = [];
    
    linearCoefficients = regress(edgePoints', [ones(length(times),1) times']); % fit a line
    invasionSpeed(i) = linearCoefficients(2);
    
%     plot(times,edgePoints)
%     hold on
%     plot(times, linearCoefficients(1) + linearCoefficients(2)*times)
%     hold off

    v = population.trait_variance;
    maxVariance(i) = v(ceil(length(v)/2), end); % maximum of v occurs at the center
end

figure,
plot(gradient, invasionSpeed, 'Color', purple, 'LineWidth', 1.5)
xlabel('Gradient of the Optimal Trait  $[\mathtt{Q} / \mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel('Invasion Wave Speed $[\mathtt{X}/\mathtt{T}]$','Interpreter','latex','FontSize', 12);
grid on


figure,
plot(gradient, amplitude, 'Color', purple, 'LineWidth', 1.5)
xlabel('Optimal Trait Gradient $[\mathtt{Q} / \mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel('Invasion Wave Amplitude $[\mathtt{N}/\mathtt{X}]$','Interpreter','latex','FontSize', 12);
grid on



%% Checking the steady state amplitude of n an v with respect to the gradient 
figure,
plot(gradient, maxVariance)
hold on
maxVarianceVsGradient = regress(maxVariance,  [ones(length(gradient),1) gradient]); % fit a line
plot(gradient, maxVarianceVsGradient(1) + maxVarianceVsGradient(2)*gradient)
hold off

v_star = 0.1:0.01:20;
V_s = modelParameters.V_s;
K = modelParameters.K;
R = modelParameters.R;
D = modelParameters.D;
V_u = modelParameters.V_u;
n_star = K .* sqrt( (v_star + V_u)./ V_u ) .* ( 1 - v_star ./ (2*V_s*R) );
figure
plot((v_star-maxVarianceVsGradient(1))./ maxVarianceVsGradient(2) , n_star)

v_star = zeros(1,length(gradient));
for i = 1 : length(gradient)
    f = @(v) ( 5/V_s * v^3 + (4/V_s*V_u - 2*R(1)) * v^2 - 8*D(1)*gradient(i)^2 * v - 8*D(1)*V_u*gradient(i)^2 );
    v_star(i) = fsolve(f, maxVariance(i));
end 
figure,
plot(gradient, v_star);
