clear; clc; close all;

% =========================================================================
% SECTION 1: Load and Visualize Data 
% =========================================================================
load('exp_growth.mat');

% The variables are 'time' and 'Npop'
disp('Variables loaded from exp_growth.mat:');
whos

figure;
scatter(time, Npop, 'filled', ...
        'MarkerFaceAlpha', 0.5, ...
        'MarkerEdgeColor', 'b');
xlabel('Time');
ylabel('Population Size');
title('Observed Data and Non-Linear Fit');
grid on;
hold off;
saveas(gcf, 'scatter_expgrowth.png');

model = @(params, t) params(1) * exp(params(2) * t);
initial_params = [1, 0.1]; % Initial guess for [N0, r]

% non-linear least squares fitting reference: https://www.mathworks.com/help/stats/nlinfit.html
[beta0, R, J, CovB, MSE] = nlinfit(time, Npop, model, initial_params);

% Extract best-fit parameters
N0_fit = beta0(1);
r_fit = beta0(2);
% Generate fitted data
Npop_fit = model(beta0, time);

figure;
scatter(time, Npop, 'filled', ...
            'MarkerFaceAlpha', 0.5, ...
            'MarkerEdgeColor', 'b');
hold on;
plot(time, Npop_fit, 'r-', 'LineWidth', 2);
xlabel('Time');
ylabel('Population Size');
title('Exponential Growth Model Fit');
legend('Observed Data', 'Fitted Curve');
grid on;
hold off;
saveas(gcf, 'fit_expgrowth.png');

% Calculate residuals
residuals = Npop - Npop_fit;

% Plot residuals
figure;
histogram(residuals, 'BinMethod', 'sturges', 'FaceColor', 'b',...
                 'EdgeColor', 'w', 'Normalization', 'pdf');
hold on;
x_vals = linspace(min(residuals)-2, max(residuals)+2, 100);
pd = fitdist(residuals, 'Normal');
y_vals = pdf(pd, x_vals);
plot(x_vals, y_vals, 'r--', 'LineWidth', 0.5);
xlabel('Residuals');
ylabel('Probability Density');
legend('Normal Distribution Fit');
grid on;
hold off;
saveas(gcf, 'residuals_hist_expgrowth.png');

% Analyze residuals
mean_res = mean(residuals);
std_res = std(residuals);
fprintf('Mean of residuals: %.4f\n', mean_res);
fprintf('Standard deviation of residuals: %.4f\n', std_res);

% Display best-fit parameters
fprintf('Best-fit parameters:\n');
fprintf('N0 = %.4f\n', N0_fit);
fprintf('r = %.4f\n', r_fit);

% Calculate 95% confidence intervals for parameters, using nlpredci and nlparci 
% reference: https://www.mathworks.com/help/stats/nlpredci.html
% and https://www.mathworks.com/help/stats/nlparci.html
[ypred, delta] = nlpredci(model, time, beta0, residuals,'Covar', CovB);
ci_N0 = nlparci(beta0, R, 'Covar', CovB);
ci_r = nlparci(beta0, R, 'Covar', CovB);

fprintf('95%% Confidence Intervals:\n');
fprintf('N0: [%.4f, %.4f]\n', ci_N0(1,1), ci_N0(1,2));
fprintf('r: [%.4f, %.4f]\n', ci_r(2,1), ci_r(2,2));

sup_lim = ypred + delta;
inf_lim = ypred - delta;

figure;
scatter(time, Npop, 'filled', ...
            'MarkerFaceAlpha', 0.5, ...
            'MarkerEdgeColor', 'b');
hold on;
plot(time, Npop_fit, 'r-', 'LineWidth', 2);
plot(time, [sup_lim, inf_lim], 'k--', 'LineWidth', 1);
xlabel('Time');
ylabel('Population Size');
title('Exponential Growth Model Fit with 95% Confidence Intervals');
legend('Observed Data', 'Fitted Curve', '95% Confidence Intervals');
grid on;
hold off;
saveas(gcf, 'fit_expgrowth_CI.png');





