clear; clc; close all;

% =========================================================================
% SECTION 1: Define Data and Model
% =========================================================================
% Data vectors 
x = [1, 2, 3, 5, 7, 10]';
y = [109, 149, 149, 191, 213, 224]';
ey = [10, 10, 2, 2, 2, 2]'; % Measurement errors for y

% Define the exponential function to be fitted
model = @(b,x) b(1) * (1 - exp(-b(2) * x));

% Provide an initial guess for the parameters [b1, b2].
b0 = [150, 0.1];

% =========================================================================
% SECTION 2: Perform Unweighted and Weighted Fits
% =========================================================================

% --- Standard (Unweighted) Non-Linear Fit ---
% All data points are treated with equal importance.
mdl_unweighted = fitnlm(x, y, model, b0);

% --- Weighted Non-Linear Fit ---
w = 1./(ey.^2);
mdl_weighted = fitnlm(x, y, model, b0, 'Weights', w);

% =========================================================================
% SECTION 3: Display and Compare Results
% =========================================================================

% Display the statistical summary for both models
disp('--- 1. Standard (Unweighted) Fit Results ---');
disp(mdl_unweighted);

disp('--- 2. Weighted Fit Results ---');
disp(mdl_weighted);

% --- Create a comparison plot ---
figure;
hold on;
box on;
grid on;

% Plot original data with error bars
errorbar(x, y, ey, 'bo', 'LineWidth', 1, 'DisplayName', 'Observed Data');

% Generate a smooth curve for plotting the models
x_curve = linspace(0, 12, 250);
y_unweighted = predict(mdl_unweighted, x_curve');
y_weighted = predict(mdl_weighted, x_curve');

% Plot the two fitted curves
plot(x_curve, y_unweighted, 'r-', 'LineWidth', 2, 'DisplayName', 'Unweighted Fit');
plot(x_curve, y_weighted, 'g--', 'LineWidth', 2, 'DisplayName', 'Weighted Fit');

% Final plot formatting
title('Comparison of Unweighted and Weighted Non-Linear Regression');
xlabel('x');
ylabel('y');
legend('Location', 'southeast');
ylim([80 240]);
hold off;
saveas(gcf, 'comparison_weighted_unweighted.png');