clear, close all, clc

%%=============================================================================
%  PROCESS EASTERN U.S. DATA
%==============================================================================

% Load and Sort Data
disp('--- Processing Eastern U.S. Data ---');
east_data = readtable('sediment_eastUS.txt');
east_data = sortrows(east_data, 'age'); % Sort by age

% Create a Scatter Plot

figure;
scatter(east_data.age, east_data.depth, 'filled', ...
            'MarkerFaceAlpha', 0.5, ...
            'MarkerEdgeColor', 'b');
xlabel('Age (kyr)');
ylabel('Depth (m)');
title('Sediment Age vs. Depth (Eastern U.S.)');
grid on;
saveas(gcf, 'scatter_eastUS.png');

% Compute Correlation Coefficient
[R, P] = corrcoef(east_data.age, east_data.depth);
r_value = R(1, 2);
p_value = P(1, 2);
fprintf('Pearson correlation coefficient (r): %.4f\n', r_value);
fprintf('P-value: %.4f\n', p_value);
if p_value < 0.05
    disp('The correlation is statistically significant (p < 0.05).');
else
    disp('The correlation is not statistically significant (p >= 0.05).');
end

% Fit a Linear Model and Plot with Confidence Bounds
[p, S] = polyfit(east_data.age, east_data.depth, 1);
depth_fit = polyval(p, east_data.age);
[depth_pred, delta] = polyconf(p, east_data.age, S);
sup_lim = depth_pred + delta;
inf_lim = depth_pred - delta;

figure;
hold on;
scatter(east_data.age, east_data.depth, 'filled', 'Color', 'b',...
         'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'b', 'DisplayName', 'Data');
plot(east_data.age, depth_fit, 'r-', 'LineWidth', 2,...
            'DisplayName', 'Best-Fit Line');
plot(east_data.age, [sup_lim, inf_lim], 'r--', 'LineWidth', 1,...
        'DisplayName', '95% Confidence Bounds');
% fill the area between confidence bounds along the whole fit line
fill([east_data.age; flipud(east_data.age)], [sup_lim; flipud(inf_lim)], 'r',...
             'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
xlabel('Age (kyr)');
ylabel('Depth (m)');
legend('Location', 'northwest');
title('Linear Fit with 95% Confidence Bounds (Eastern U.S.)');
grid on;
hold off;
fprintf('The linear model is: depth = (%.4f) * age + (%.4f)\n', p(1), p(2));
saveas(gcf, 'linear_fit_eastUS.png');

% Plot Histogram of Residuals (Eastern U.S.)
residuals_east = east_data.depth - depth_fit;
figure;
% sturges method for bin width selection and pdf normalization for better visualization 
% reference: https://www.mathworks.com/help/matlab/ref/matlab.graphics.chart.primitive.histogram.html
histogram(residuals_east, 'BinMethod', 'sturges', 'Normalization', 'pdf',...
            'FaceColor', 'b', 'EdgeColor', 'w', 'HandleVisibility', 'off');
hold on;
x_vals = linspace(min(residuals_east)-2, max(residuals_east)+2, 100);
pd = fitdist(residuals_east, 'Normal');
y_vals = pdf(pd, x_vals);
plot(x_vals, y_vals, 'r--', 'LineWidth', 0.5);
xlabel('Residuals (m)');
ylabel('Probability Density');
title('Histogram of Residuals (Eastern U.S.)');
legend('Normal Distribution Fit');
grid on;
hold off;
saveas(gcf, 'residuals_hist_eastUS.png');

% Test Residuals for Normal Distribution (Eastern U.S.)
[h, p_kstest] = kstest(residuals_east);
if h == 0
    disp('The residuals appear to be normally distributed (K-S test, p > 0.05).');
else
    disp('The residuals do not appear to be normally distributed (K-S test, p <= 0.05).');
end
fprintf('K-S test p-value: %.4f\n', p_kstest);

% Compute Coefficient of Determination (R^2) for Eastern U.S. Data
TSS_east = sum((east_data.depth - mean(east_data.depth)).^2);
RSS_east = sum(residuals_east.^2);
R2_east = 1 - (RSS_east / TSS_east);
fprintf('Coefficient of Determination (R^2) for Eastern U.S.: %.4f\n\n', R2_east);


%%=============================================================================
%  PROCESS WESTERN U.S. DATA
%==============================================================================

% Load and Sort Data
disp('--- Processing Western U.S. Data ---');
west_data = readtable('sediment_westUS.txt');
west_data = sortrows(west_data, 'age');

% Create a Scatter Plot (Western U.S.)
figure;
scatter(west_data.age, west_data.depth, 'filled',...
        'MarkerFaceColor', 'green', 'MarkerFaceAlpha', 0.5,...
         'MarkerEdgeColor', 'green');
xlabel('Age (kyr)');
ylabel('Depth (m)');
title('Sediment Age vs. Depth (Western U.S.)');
grid on;
saveas(gcf, 'scatter_westUS.png');

% Compute Correlation Coefficient (Western U.S.)
[R_west, P_west] = corrcoef(west_data.age, west_data.depth);
r_value_west = R_west(1, 2);
p_value_west = P_west(1, 2);
fprintf('Pearson correlation coefficient (r): %.4f\n', r_value_west);
fprintf('P-value: %.4f\n', p_value_west);
if p_value_west < 0.05
    disp('The correlation is statistically significant (p < 0.05).');
else
    disp('The correlation is not statistically significant (p >= 0.05).');
end

% Fit a Linear Model and Plot with Confidence Bounds (Western U.S.)
[p_west, S_west] = polyfit(west_data.age, west_data.depth, 1);
depth_fit_west = polyval(p_west, west_data.age);
[depth_pred_west, delta_west] = polyconf(p_west, west_data.age, S_west);
sup_lim_west = depth_fit_west + delta_west;
inf_lim_west = depth_fit_west - delta_west;



figure;
hold on;
scatter(west_data.age, west_data.depth, 'filled', 'Color', 'green',...
     'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'green', 'DisplayName', 'Data');
plot(west_data.age, depth_fit_west, 'r-', 'LineWidth', 2,...
        'DisplayName', 'Best-Fit Line');
plot(west_data.age, [sup_lim_west, inf_lim_west], 'r--', 'LineWidth', 1,...
         'DisplayName', '95% Confidence Bounds');
fill([west_data.age; flipud(west_data.age)], [sup_lim_west; flipud(inf_lim_west)], 'r',...
             'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
xlabel('Age (kyr)');
ylabel('Depth (m)');
legend('Location', 'northwest');
title('Linear Fit with 95% Confidence Bounds (Western U.S.)');
grid on;
hold off;
fprintf('The linear model is: depth = (%.4f) * age + (%.4f)\n', p_west(1), p_west(2));
saveas(gcf, 'linear_fit_westUS.png');


% Plot Histogram of Residuals (Western U.S.)
residuals_west = west_data.depth - depth_fit_west;
figure;
histogram(residuals_west, 'BinMethod', 'sturges', 'Normalization', 'pdf',...
             'FaceColor', 'green', 'EdgeColor', 'w', 'HandleVisibility', 'off');
hold on;
x_vals_west = linspace(min(residuals_west)-2, max(residuals_west)+2, 100);
pd_west = fitdist(residuals_west, 'Normal');
y_vals_west = pdf(pd_west, x_vals_west);
plot(x_vals_west, y_vals_west, 'r--', 'LineWidth', 0.5);
xlabel('Residuals (m)');
ylabel('Probability Density');
title('Histogram of Residuals (Western U.S.)');
legend('Normal Distribution Fit');
grid on;
hold off;
saveas(gcf, 'residuals_hist_westUS.png');

% Test Residuals for Normal Distribution (Western U.S.)
[h_west, p_kstest_west] = kstest(residuals_west);
if h_west == 0
    disp('The residuals appear to be normally distributed (K-S test, p > 0.05).');
else
    disp('The residuals do not appear to be normally distributed (K-S test, p <= 0.05).');
end
fprintf('K-S test p-value: %.4f\n', p_kstest_west);

% Compute Coefficient of Determination (R^2) for Western U.S. Data
TSS_west = sum((west_data.depth - mean(west_data.depth)).^2);
RSS_west = sum(residuals_west.^2);
R2_west = 1 - (RSS_west / TSS_west);
fprintf('Coefficient of Determination (R^2) for Western U.S.: %.4f\n\n', R2_west);

%%=============================================================================
%  BOOTSTRAP ANALYSIS FOR EASTERN U.S. DATA
%==============================================================================

disp('--- Bootstrap Analysis for Eastern U.S. Data ---');
% Define number of bootstrap samples
N_bootstrap = 1000;

disp('--- Performing Bootstrap Analysis for Pearson r ---');
N_boot = 1000;

% Define the function to compute correlation for bootstrapping
corr_func = @(age, depth) corr(age, depth);

% Bootstrap for Eastern US data
r_boot_east = bootstrp(N_boot, corr_func, east_data.age, east_data.depth);

% Bootstrap for Western US data
r_boot_west = bootstrp(N_boot, corr_func, west_data.age, west_data.depth);

% Plot the histograms
figure;
subplot(1, 2, 1);
histogram(r_boot_east, 'FaceColor', 'blue', 'EdgeColor', 'w');
title('Bootstrap of r-values (Eastern US)');
xlabel('Pearson r');
ylabel('Frequency');
grid on;

subplot(1, 2, 2);
histogram(r_boot_west, 'FaceColor', 'green', 'EdgeColor', 'w');
title('Bootstrap of r-values (Western US)');
xlabel('Pearson r');
ylabel('Frequency');
grid on;
saveas(gcf, 'bootstrap_rvalues.png');

% Display insights
mean_r_east = mean(r_boot_east);
std_r_east = std(r_boot_east);
mean_r_west = mean(r_boot_west);
std_r_west = std(r_boot_west);

fprintf('Eastern US: Mean r = %.4f, Std Dev = %.4f\n', mean_r_east, std_r_east);
fprintf('Western US: Mean r = %.4f, Std Dev = %.4f\n', mean_r_west, std_r_west);

disp('--- Bootstrapping Regression Coefficients ---');

% Define the polyfit function for bootstrapping
polyfit_func = @(age, depth) polyfit(age, depth, 1);

% Bootstrap for Eastern US data
coeffs_boot_east = bootstrp(N_boot, polyfit_func, east_data.age, east_data.depth);
slopes_east = coeffs_boot_east(:, 1);
intercepts_east = coeffs_boot_east(:, 2);

% Bootstrap for Western US data
coeffs_boot_west = bootstrp(N_boot, polyfit_func, west_data.age, west_data.depth);
slopes_west = coeffs_boot_west(:, 1);
intercepts_west = coeffs_boot_west(:, 2);

% Calculate 95% confidence intervals
ci_slope_east = prctile(slopes_east, [2.5, 97.5]);
ci_intercept_east = prctile(intercepts_east, [2.5, 97.5]);
ci_slope_west = prctile(slopes_west, [2.5, 97.5]);
ci_intercept_west = prctile(intercepts_west, [2.5, 97.5]);

fprintf('Eastern US 95%% CI for Slope: [%.4f, %.4f]\n',...
             ci_slope_east(1), ci_slope_east(2));
fprintf('Eastern US 95%% CI for Intercept: [%.4f, %.4f]\n',...
             ci_intercept_east(1), ci_intercept_east(2));
fprintf('Western US 95%% CI for Slope: [%.4f, %.4f]\n',...
             ci_slope_west(1), ci_slope_west(2));
fprintf('Western US 95%% CI for Intercept: [%.4f, %.4f]\n\n',...
             ci_intercept_west(1), ci_intercept_west(2));

% Plot the distributions
figure;
subplot(2, 2, 1);
histogram(slopes_east, 'FaceColor', 'blue', 'EdgeColor', 'w', 'HandleVisibility', 'off');
xline(ci_slope_east(1), 'r--', 'LineWidth', 1.5,...
            'label', sprintf('%.4f', ci_slope_east(1)),...
            'LabelVerticalAlignment', 'middle',...
            'LabelHorizontalAlignment', 'right');
xline(ci_slope_east(2), 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off',...
                     'label', sprintf('%.4f', ci_slope_east(2)),...
                     'LabelVerticalAlignment', 'middle',...
                     'LabelHorizontalAlignment', 'left');
legend('95% CI');
title('Slope Distribution (East US)');
xlabel('Slope'); ylabel('Frequency'); grid on;

subplot(2, 2, 2);
histogram(intercepts_east, 'FaceColor', 'blue', 'EdgeColor', 'w', 'HandleVisibility', 'off');
xline(ci_intercept_east(1), 'r--', 'LineWidth', 1.5, ...
      'label', sprintf('%.4f', ci_intercept_east(1)),...
      'LabelVerticalAlignment', 'middle',...
      'LabelHorizontalAlignment', 'right');
xline(ci_intercept_east(2), 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off', ...
      'label', sprintf('%.4f', ci_intercept_east(2)),...
      'LabelVerticalAlignment', 'middle',...
      'LabelHorizontalAlignment', 'left');
legend('95% CI');
title('Intercept Distribution (East US)');
xlabel('Intercept'); ylabel('Frequency'); grid on;

subplot(2, 2, 3);
histogram(slopes_west, 'FaceColor', 'green', 'EdgeColor', 'w', 'HandleVisibility', 'off');
xline(ci_slope_west(1), 'r--', 'LineWidth', 1.5,...
        'label', sprintf('%.4f', ci_slope_west(1)),...
        'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'right');
xline(ci_slope_west(2), 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off', ...
        'label', sprintf('%.4f', ci_slope_west(2)),...
        'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'left');
legend('95% CI');
title('Slope Distribution (West US)');
xlabel('Slope'); ylabel('Frequency'); grid on;

subplot(2, 2, 4);
histogram(intercepts_west, 'FaceColor', 'green', 'EdgeColor', 'w', 'HandleVisibility', 'off');
xline(ci_intercept_west(1), 'r--', 'LineWidth', 1.5,...
        'label', sprintf('%.4f', ci_intercept_west(1)),...
        'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'right');
xline(ci_intercept_west(2), 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off', ...
        'label', sprintf('%.4f', ci_intercept_west(2)),...
        'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'left');
legend('95% CI');
title('Intercept Distribution (West US)');
xlabel('Intercept'); ylabel('Frequency'); grid on;
saveas(gcf, 'bootstrap_coefficients.png');


disp('--- Removing Outliers from Western U.S. Data and Refitting ---');

% Remove outlier(s) based on the 'depth' variable
[~, outlier_indices] = rmoutliers(west_data.depth);
west_data_clean = west_data(~outlier_indices, :);

fprintf('Removed %d outlier(s).\n', sum(outlier_indices));

% Re-compute Correlation
[R_clean, P_clean] = corrcoef(west_data_clean.age, west_data_clean.depth);
fprintf('New Pearson r: %.4f\n', R_clean(1,2));
fprintf('New P-value: %.4f\n', P_clean(1,2));

% Re-compute R^2
residuals_clean = west_data_clean.depth - polyval(polyfit(west_data_clean.age, west_data_clean.depth, 1), west_data_clean.age);
TSS_clean = sum((west_data_clean.depth - mean(west_data_clean.depth)).^2);
RSS_clean = sum(residuals_clean.^2);
R2_clean = 1 - (RSS_clean / TSS_clean);
fprintf('New R^2: %.4f\n', R2_clean);

% Re-fit and plot the model
p_clean = polyfit(west_data_clean.age, west_data_clean.depth, 1);
depth_fit_clean = polyval(p_clean, west_data_clean.age);
[depth_pred_clean, delta_clean] = polyconf(p_clean, west_data_clean.age, S_west);
sup_lim_clean = depth_pred_clean + delta_clean;
inf_lim_clean = depth_pred_clean - delta_clean;

figure;
hold on;
scatter(west_data_clean.age, west_data_clean.depth, 'filled',...
        'Color', 'green', 'MarkerFaceAlpha', 0.6,...
        'MarkerEdgeColor', 'green', 'DisplayName', 'Data');
plot(west_data_clean.age, depth_fit_clean, 'r-', 'LineWidth', 2,...
        'DisplayName', 'Best-Fit Line');
plot(west_data_clean.age, [sup_lim_clean, inf_lim_clean], 'r--',...
         'LineWidth', 1, 'DisplayName', '95% Confidence Bounds');

% fill the area between confidence bounds along the whole fit line
fill([west_data_clean.age; flipud(west_data_clean.age)],...
         [sup_lim_clean; flipud(inf_lim_clean)], 'r', ...
         'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
xlabel('Age (kyr)');
ylabel('Depth (m)');
legend('Location', 'northwest');
title('Linear Fit with 95% Confidence Bounds (Western U.S.)');
grid on;
hold off;
saveas(gcf, 'linear_fit_westUS_clean.png');
fprintf('The linear model is: depth = (%.4f) * age + (%.4f)\n', p_clean(1), p_clean(2));