%%
%
%       CERI 8104/CIVL 8126 Data Analysis in Geophysics
%
%       fit well draw down - invert for S and T
%       Solution: S = 2*1e-4; [] unitless
%                     T = 14000*ft2/day;
%       from Lohman 1975, USGS open file report
% 
%       @author thgoebel, CERI U of Memphis
%       Modified by A Martinez-Coto, Fall 2025
%
%
clear; clc; close all;

ft2m      = 0.3048;
fontsize = 16;
dy2sc    = 3600*24;

%%=======================1===============================
%             set params and files
%=======================================================
file_in = 'draw_down_data.csv';
% known params: switch here between 200, 400 and 800 feet
dist = 200;
r        = dist*ft2m;
Q       = 96000*ft2m^3/dy2sc; %ft3/day -> m3/s
sprintf( "pump rate : %.3f (m3/s)", Q)

% initial guess
S0 = 1e-4; T0 = 1e-3; %m2/s
par_lb  = [1e-5, 0.001];
par_ub  = [5*1e-3, 0.1];

%%=======================2===============================
%          load data  and convert to SI units
%=======================================================
data = readtable( file_in, 'HeaderLines', 2);
data.Time = data.Time*60; %/60/24;%*60
if r == 200*ft2m
    data.drawdown= data.drawdown200*ft2m;
elseif r == 400*ft2m
    data.drawdown= data.drawdown400*ft2m;
elseif r == 800*ft2m
    data.drawdown= data.drawdown800*ft2m;
end

figure('Name', 'Data Visualization');
subplot(1,3,1) 
scatter( data.Time, data.drawdown, 'bo', 'filled', ...
         'MarkerFaceAlpha', 0.5, ...
         'MarkerEdgeColor', 'b')
xlabel('Time (s)', 'FontSize', fontsize)
ylabel('Water Level (m)', 'FontSize', fontsize)
grid on
subplot(1,3,2)
scatter( data.Time, data.drawdown, 'ro', 'filled', ...
         'MarkerFaceAlpha', 0.5, ...
         'MarkerEdgeColor', 'r')
xlabel('Time (s)', 'FontSize', fontsize)
ylabel('Water Level (m)', 'FontSize', fontsize)
% I decided to use set instead of semilogx or semilogy because it gives more 
% flexibility when setting the markers properties and colors
set(gca, 'XScale', 'log', 'YScale', 'log');
grid on
subplot(1,3,3)
scatter( data.Time, data.drawdown, 'go', 'filled', ...
         'MarkerFaceAlpha', 0.5, ...
         'MarkerEdgeColor', 'g')
set(gca, 'XScale', 'log'); 
xlabel('Time (s)', 'FontSize', fontsize)
ylabel('Water Level (m)', 'FontSize', fontsize)
grid on
sgtitle(sprintf('Water Level Drawdown at %.2f m', r), 'FontSize', fontsize, 'FontWeight', 'bold');

% %%=======================3===============================
% %                      invert for S and T
% %=======================================================
modelFct = @(par,t)well_fct( t, r, Q, par);
par_hat = nlinfit( data.Time, data.drawdown, modelFct, [S0, T0]);

sprintf( "Unconstrainted LSQ: S=%d, T=%.2f", par_hat)
sprintf( "Hydrau. Diffusivity: %.2f (m2/s)", (par_hat(2)/par_hat(1)))

% prediction with best-fit model
a_h_hat = well_fct(  data.Time, r, Q, par_hat);

% constrained lsqcurvefit
par_hat_constr    = lsqcurvefit(modelFct, [S0, T0],...
                                 data.Time, data.drawdown,...
                                 par_lb, par_ub);
sprintf( "Constrainted LSQ: S=%d, T=%.2f", par_hat_constr)
% prediction with best-fit model
a_h_hat_constr = well_fct(  data.Time, r, Q, par_hat_constr);

% %=======================4===============================
% %         plot data and best LSQ solution
% %=======================================================
% make a two panel figure
figure('Name', 'Data and Theis Model Fit');
subplot(1,2,1)
scatter( data.Time, data.drawdown, 'bo', 'filled', ...
         'MarkerFaceAlpha', 0.5, ...
         'MarkerEdgeColor', 'b')
hold on
plot( data.Time, a_h_hat, 'r-', 'LineWidth', 2)
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Time (s)', 'FontSize', fontsize)
ylabel('Water Level (m)', 'FontSize', fontsize)
title(sprintf('S=%.8e, T=%.8f', par_hat(1), par_hat(2)), 'FontSize', fontsize)
legend('Data', 'Unconstrained Fit', 'FontSize', fontsize)
grid on
hold off
subplot(1,2,2)
scatter( data.Time, data.drawdown, 'bo', 'filled', ...
         'MarkerFaceAlpha', 0.5, ...
            'MarkerEdgeColor', 'b')
hold on
plot( data.Time, a_h_hat_constr, 'g-', 'LineWidth', 2)
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Time (s)', 'FontSize', fontsize)
ylabel('Water Level (m)', 'FontSize', fontsize)
title(sprintf('T=%.8e, S=%.8f', par_hat_constr(1), par_hat_constr(2)), 'FontSize', fontsize)
legend('Data', 'Constrained Fit', 'FontSize', fontsize)
grid on
hold off
sgtitle(sprintf('Water Level Drawdown at %.2f m', r), 'FontSize', fontsize, 'FontWeight', 'bold');

%=======================5===============================
%         plot residuals and do stats
%=======================================================
% Analyze residuals
res_unconstr = data.drawdown - a_h_hat;
res_constr = data.drawdown - a_h_hat_constr;

mean_res_unconstr = mean(res_unconstr);
std_res_unconstr = std(res_unconstr);
fprintf('Mean of unconstrained residuals: %.4f\n', mean_res_unconstr);
fprintf('Standard deviation of unconstrained residuals: %.4f\n', std_res_unconstr);

mean_res_constr = mean(res_constr);
std_res_constr = std(res_constr);
fprintf('Mean of constrained residuals: %.4f\n', mean_res_constr);
fprintf('Standard deviation of constrained residuals: %.4f\n', std_res_constr);

% plot residuals to see if there are any outliers
figure('Name', 'Residual Analysis');
subplot(2,1,1);
binsize = 0.01;
numBins = ceil((max(res_unconstr)-min(res_unconstr))/binsize);

histogram(res_unconstr, 'FaceColor', 'b', 'EdgeColor', 'w',...
            'Normalization', 'pdf', 'NumBins', numBins);
hold on;
x_vals = linspace(min(res_unconstr)-0.5, max(res_unconstr)+0.1, 100);
pd = fitdist(res_unconstr, 'Normal');
y_vals = pdf(pd, x_vals);
plot(x_vals, y_vals, 'r--', 'LineWidth', 0.5);
xline(mean_res_unconstr + 4*std_res_unconstr, 'k--', 'LineWidth', 0.5, 'Label', '4\sigma',...
     'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'right');
xline(mean_res_unconstr - 4*std_res_unconstr, 'k--', 'LineWidth', 0.5, 'Label', '-4\sigma',...
    'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'left');
xlabel('Residuals');
ylabel('Probability Density');
legend('Residuals Unconstrained', 'Normal Distribution Fit');
grid on;
subplot(2,1,2);
histogram(res_constr, 'FaceColor', 'g', 'EdgeColor', 'w',...
            'Normalization', 'pdf', 'NumBins', numBins);
hold on;
x_vals_constr = linspace(min(res_constr)-0.5, max(res_constr)+0.1, 100);
pd_constr = fitdist(res_constr, 'Normal');
y_vals_constr = pdf(pd_constr, x_vals_constr);
plot(x_vals_constr, y_vals_constr, 'r--', 'LineWidth', 0.5);
xline(mean_res_constr + 4*std_res_constr, 'k--', 'LineWidth', 0.5, 'Label', '4\sigma',...
     'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'right');
xline(mean_res_constr - 4*std_res_constr, 'k--', 'LineWidth', 0.5, 'Label', '-4\sigma',...
    'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'left');
xlabel('Residuals');
ylabel('Probability Density');
legend('Residuals Histogram', 'Normal Distribution Fit');
sgtitle('Residual Analysis', 'FontSize', fontsize, 'FontWeight', 'bold');
grid on;
hold off;

%=======================6===============================
%         Bootstrap to estimate parameter uncertainty
%=======================================================

% from here I am using the constrained parameters only

nBoot = 1000; % Number of bootstrap samples
warning('off');
coeff_boot = bootstrp(nBoot, @(x,y) lsqcurvefit(modelFct, [S0, T0],...
                         x, y, par_lb, par_ub), data.Time, data.drawdown);
ci_S = prctile(coeff_boot(:,1), [2.5 97.5]);
ci_T = prctile(coeff_boot(:,2), [2.5 97.5]);

% Plot histogram of coefficients from bootstrapping
figure('Name', 'Bootstrap Parameter Distributions');
subplot(1,2,1);
histogram(coeff_boot(:,1), 'FaceColor', 'b', 'EdgeColor', 'w');
xline(ci_S(1), 'r--', 'LineWidth', 1.5, 'Label', sprintf('%.8f', ci_S(1)),...
        'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'right');
xline(ci_S(2), 'r--', 'LineWidth', 1.5, 'Label', sprintf('%.8f', ci_S(2)),...
        'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'left');
xlabel('Storage Coefficient S', 'FontSize', fontsize);
ylabel('Frequency', 'FontSize', fontsize);
title('Bootstrap Distribution of S', 'FontSize', fontsize);
grid on;
subplot(1,2,2);
histogram(coeff_boot(:,2), 'FaceColor', 'r', 'EdgeColor', 'w');
xline(ci_T(1), 'b--', 'LineWidth', 1.5, 'Label', sprintf('%.6f', ci_T(1)),...
        'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'right');
xline(ci_T(2), 'b--', 'LineWidth', 1.5, 'Label', sprintf('%.6f', ci_T(2)),...
       'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'left');
xlabel('Transmissivity T', 'FontSize', fontsize);
ylabel('Frequency', 'FontSize', fontsize);
title('Bootstrap Distribution of T', 'FontSize', fontsize);
grid on;

%======================================================================= 
%           remove outliers based on the highest residuals
%=======================================================================

% Define a threshold for outlier, in this case 4 standard deviations 
threshold = 4 * std_res_constr;

% Identify outliers
outlier_indices = abs(res_constr) > threshold;

if any(outlier_indices)
    fprintf('Outliers detected at indices: %s\n', mat2str(find(outlier_indices)));
else
    fprintf('No outliers detected based on the threshold of %.4f\n', threshold);
end

% Remove outliers from data
data_cleaned = data(~outlier_indices, :);
% Refit the model to the cleaned data
par_hat_cleaned = lsqcurvefit(modelFct, [S0, T0],...
                                 data_cleaned.Time, data_cleaned.drawdown,...
                                 par_lb, par_ub);
sprintf( "After Outlier Removal - Constrainted LSQ: S=%d, T=%.2f", par_hat_cleaned)
sprintf( "Hydrau. Diffusivity: %.2 (m2/s)", (par_hat_cleaned(2)/par_hat_cleaned(1)))
% prediction with best-fit model
a_h_hat_cleaned = well_fct(  data_cleaned.Time, r, Q, par_hat_cleaned);
% Plot the cleaned data and new fit
figure('Name', 'Data and Theis Model Fit (After Outlier Removal)');
scatter( data_cleaned.Time, data_cleaned.drawdown, 'bo', 'filled', ...
         'MarkerFaceAlpha', 0.5, ...
         'MarkerEdgeColor', 'b')
hold on
plot( data_cleaned.Time, a_h_hat_cleaned, 'r-', 'LineWidth', 2)
xlabel('Time (s)', 'FontSize', fontsize)
ylabel('Water Level (m)', 'FontSize', fontsize)
set(gca, 'XScale', 'log', 'YScale', 'log');
title(sprintf('Theis Model Fit at r=%.2f m (Outlier Removal)', r), 'FontSize', fontsize)
legend('Cleaned Data', 'Constrained Fit', 'FontSize', fontsize)
grid on
hold off

residuals_cleaned = data_cleaned.drawdown - a_h_hat_cleaned;
mean_res_cleaned = mean(residuals_cleaned);
std_res_cleaned = std(residuals_cleaned);

%==============================================================================
% Display results summary
%==============================================================================
fprintf('\n=========== RESULTS SUMMARY at r=%.2f m ====================\n', r);
fprintf('Initial Fit (With Outliers):\n');
fprintf('  Storage Coefficient Unconstrained S: %.10f []\n', par_hat(1));
fprintf('  Transmissivity Unconstrained T: %.10f (m2/s)\n', par_hat(2));
fprintf('  Hydraulic Diffusivity Unconstrained D: %.10f (m2/s)\n', (par_hat(2)/par_hat(1))); 
fprintf('  Storage Coefficient Constrained S: %.10f []\n', par_hat_constr(1));
fprintf('  Transmissivity Constrained T: %.10f (m2/s)\n', par_hat_constr(2));
fprintf('  Hydraulic Diffusivity Constrained D: %.10f (m2/s)\n', (par_hat_constr(2)/par_hat_constr(1)));
fprintf('  Mean of Residuals (Unconstrained): %.10f\n', mean_res_unconstr);
fprintf('  Standard Deviation of Residuals (Unconstrained): %.10f\n', std_res_unconstr);
fprintf('  Mean of Residuals (Constrained): %.10f\n', mean_res_constr);
fprintf('  Standard Deviation of Residuals (Constrained): %.10f\n', std_res_constr);
fprintf('\nAfter Outlier Removal:\n');
fprintf('  Storage Coefficient S: %.10f []\n', par_hat_cleaned(1));
fprintf('  Transmissivity T: %.10f (m2/s)\n', par_hat_cleaned(2));
fprintf('  Hydraulic Diffusivity D: %.10f (m2/s)\n', (par_hat_cleaned(2)/par_hat_cleaned(1)));
fprintf('  Mean of Residuals: %.10f\n', mean_res_cleaned);
fprintf('  Standard Deviation of Residuals: %.10f\n', std_res_cleaned);
fprintf('======================================================================\n');

% %=======================================================
% %             new fct definitions
% %=======================================================

function h = well_fct(  a_t, r, Q, par)
    % takes global variables g_r and g_Q
    S = par(1);
    T = par(2);
    u = (r^2 * S) ./ (4 * T * a_t);
    h = (Q / (4 * pi * T)) * expint(u);
end
