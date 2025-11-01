%
% CERI 8104 / CIVL 8126 Data Analysis in Geophysics
%               fault roughness analysis
%      @author thgoebel, CERI U of Memphis
%
%
clear all
clf
close 

%=======================1==================================================
%                  parameters
%==========================================================================
data_dir = '/home/pequenojn/gdrive/AdoMarti/01_Projects/DataAnalysis-Fall2025/homeworks/hmwk2_faultroughness/data/'; % data directory
file_in  = "fault_roughness.txt"; % input data file

% if = 1 show pcolor of detrend surface
plot_surf = 1;
zMin = -40; % only for surface plots
zMax =  40;

%=======================2==================================================
%                   load data, gridding, detrending etc.
%==========================================================================
XYZ = load(strcat(data_dir, file_in));
%A% remove nans
XYZ =  XYZ(~isnan(XYZ(:, 3)),:);
X     =  XYZ(:, 1);
Y     =  XYZ(:, 2);
Z     =  XYZ(:, 3);

sprintf( '%s %.1f %.1f','range of topography:',min(Z), max(Z))
% spatial increments
fprintf( '%s%.2f', 'grid spacing: ', X(2) - X(1), ' m')

% interpolation to fill gaps from removed outliers
[XX, YY] = meshgrid( unique( X), unique( Y));
ZZ = griddata( X, Y, Z, XX, YY, 'linear');
%D% detrend
ZZ_detr = detrend_2d( ZZ);

% %x(mu) y(mu) z(microns)
% 0 0 475.63195

%C%visualize raw data
if plot_surf == 1
    figure(1)
    set( gcf, 'Position', [1 1 2400 800]);
    clf()
    subplot(121)
    pcolor( XX, YY, ZZ); 
    title( 'raw data')
    xlabel('X (\mu m)')
    ylabel('Y (\mu m)')
    colormap hot
    shading interp;
    cbar = colorbar('eastoutside');
    cbar.Label.String = 'Roughness';

    subplot(122)
    surf( XX, YY, ZZ); 
    title( 'raw data (3D)')
    xlabel('X (\mu m)')
    ylabel('Y (\mu m)')
    colormap hot
    shading interp;
    cbar = colorbar('eastoutside');
    cbar.Label.String = 'Roughness';
    saveas(gcf,'raw_surface.png')

    figure(2)
    clf()
    set( gcf, 'Position', [1 1 2400 800]);
    subplot(121)
    pcolor( XX, YY, ZZ_detr); 
    shading flat;
    title( 'detrended')
    xlabel('X (\mu m)')
    ylabel('Y (\mu m)')
    caxis([zMin, zMax]);
    colormap hot
    cbar = colorbar('eastoutside');
    caxis([zMin, zMax]);
    cbar.Label.String = 'Roughness';
    
    subplot(122)
    surf( XX, YY, ZZ_detr); 
    title( 'detrended (3D)')
    shading interp;
    xlabel('X (\mu m)')
    ylabel('Y (\mu m)')
    caxis([zMin, zMax]);
    colormap hot
    cbar = colorbar('eastoutside');
    caxis([zMin, zMax]);
    cbar.Label.String = 'Roughness';
    saveas(gcf,'detrended_surface.png')
end

%=======================3==================================================
%                                 stack PSDs of roughness profiles
%==========================================================================
[Ny,Nx] = size( ZZ);
dy = YY(2,1)-YY(1,1);
dx = XX(1,2)-XX(1,1);

psd_stack_x = zeros(Nx, Ny);
psd_stack_y = zeros(Ny, Nx);

%-------stack all x -profiles
for i=1:Nx-1
    profile_x = ZZ_detr(:,i);
    profile_x = profile_x - mean(profile_x); %demean
    profile_x = detrend(profile_x); 

    nfft = length(profile_x);
    fftz      = fft( profile_x, [], 1);
    psd_stack_x(i,:) = abs(fftz).^2*dx;
end

%-------stack all y -profiles
for i=1:Ny-1
    profile_y = ZZ_detr(i,:);
    profile_y = profile_y - mean(profile_y); %demean
    profile_y = detrend(profile_y);
    fftz      = fft( profile_y, [], 2);
    psd_stack_y(i,:) = abs(fftz).^2*dy;
end

% get the mean of all profiles
avg_psd_x=mean( psd_stack_x);
avg_psd_y=mean( psd_stack_y);

%=======================4==================================================
%                              plot PSD
%==========================================================================
fs = 1/dx; %sampling frequency

avg_freqx = (0:Nx-1)*(fs/Nx); %frequency vector for x-profiles
avg_freqy = (0:Ny-1)*(fs/Ny); %frequency vector for y-profiles

wavelength_x = 1./avg_freqx; %wavelength vector for x-profiles
wavelength_y = 1./avg_freqy; %wavelength vector for y-profiles

figure(3)
clf()
set( gcf, 'Position', [1 1 2400 800]);
%plot ave. power spectrum over wavelength
subplot(1,2,1)
scatter(wavelength_x, avg_psd_x, 'r', 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot(wavelength_x, avg_psd_x, 'r', 'LineWidth', 1.5);
grid on;
xlim([5, 1.5*1e3]);
xlabel('Wavelength (\mu m)');
ylabel('PSD (\mu m^3)');
set(gca, 'XScale', 'log', 'YScale', 'log');
legend('Slip Perpendicular');

subplot(1,2,2)
scatter(wavelength_y, avg_psd_y, 'b', 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot(wavelength_y, avg_psd_y, 'b', 'LineWidth', 1.5);
grid on;
xlim([5, 1.5*1e3]);
xlabel('Wavelength (\mu m)');
ylabel('PSD (\mu m^3)');
set(gca, 'XScale', 'log', 'YScale', 'log');
legend('Slip Parallel');

saveas(gcf,'psd_profiles.png');

%=======================5==================================================
%                  fit roughness exponent
%==========================================================================
fprintf( '%s\n', 'Select Roughness Exp. Fitting Range')
[x,y] = ginput(2);
fprintf( '%s%.1f%s%.1f\n', 'Selected wavelength range: ', x(1), ' to ', x(2), ' m');

% perform linear fit to log-log data that is selected by ginput
ind_x = find( wavelength_x >= min(x) & wavelength_x <= max(x));
ind_y = find( wavelength_y >= min(x) & wavelength_y <= max(x));

corner_wavelength_x = max(wavelength_x(ind_x));
corner_wavelength_y = max(wavelength_y(ind_y));

%-log-transform and lsq fit
[par_x, R_x] = polyfit( log10(wavelength_x(ind_x)), log10(avg_psd_x(ind_x)), 1);
[par_y, R_y] = polyfit( log10(wavelength_y(ind_y)), log10(avg_psd_y(ind_y)), 1);
gamma_x = par_x(1);
gamma_y = par_y(1);
Hurst_x    = (gamma_x-1)*.5;
Hurst_y    = (gamma_y-1)*.5;
sprintf(  "roughness exponent (slip perp.): %.1f, Hurst=%.2f", gamma_x, Hurst_x)
sprintf(  "roughness exponent (slip parall.): %.1f, Hurst=%.2f", gamma_y, Hurst_y)

% get fitted values for plotting
xvals = linspace(1, 2*1e3, 100);

fit_x = polyval( par_x, log10(xvals));
fit_y = polyval( par_y, log10(xvals));

figure(4)
clf()
set( gcf, 'Position', [1 1 2400 800]);
subplot(1,2,1)
scatter(wavelength_x, avg_psd_x, 'r', 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot( xvals, 10.^(fit_x), 'k--', 'LineWidth', 2);
xline(x(1), 'g--', 'LineWidth', 1.5);
xline(x(2), 'g--', 'LineWidth', 1.5);
grid on;
xlim([5, 1.5*1e3]);
xlabel('Wavelength (\mu m)');
ylabel('PSD (\mu m^3)');
set(gca, 'XScale', 'log', 'YScale', 'log');
legend('Slip Perpendicular', 'Power-law fit', 'Frequency range');

subplot(1,2,2)
scatter(wavelength_y, avg_psd_y, 'b', 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot( xvals, 10.^(fit_y), 'k--', 'LineWidth', 2);
xline(x(1), 'g--', 'LineWidth', 1.5);
xline(x(2), 'g--', 'LineWidth', 1.5);
grid on;
xlim([5, 1.5*1e3]);
xlabel('Wavelength (mu m)');
ylabel('PSD (m^3)');
set(gca, 'XScale', 'log', 'YScale', 'log');
legend('Slip Parallel', 'Power-law fit', 'Frequency range');

saveas(gcf,'psd_fits.png');

%=======================6==================================================
%                  compute rms roughness
%==========================================================================
rms_roughness = sqrt( mean( ZZ_detr(:).^2));
sprintf( '%s%.2f%s', 'RMS roughness: ', rms_roughness, ' \mu m')
sprintf( '%s%.2f%s', 'Corner wavelength: ', corner_wavelength_x)

% ===========================================================================
%                           functions
% ===========================================================================

function Z_f = detrend_2d(Z)

    %This function is written by Munther Gdeisat-The General Engineering
    %  Research Institute (GERI) at Liverpool John Moores University.
    % This program is written on 9th October 2011

    %This function is the 2D equivalent of detrend function in Matlab
    %  Z_f = DETREND(Z) removes the best plane fit trend from the
    %     data in the 2D array Z and returns the residual in the 2D array Z_f

    %Thanks for
    %    http://www.mathworks.co.uk/support/solutions/en/data/1-1AVW5/index.html?solution=1-1AVW5
    if size(Z,2) < 2
        disp('Z must be a 2D array')
        return
    end

    M = size(Z,2);
    N = size(Z,1);
    [X,Y] = meshgrid(1:M,1:N);

    %Make the 2D data as 1D vector
    Xcolv = X(:); % Make X a column vector
    Ycolv = Y(:); % Make Y a column vector
    Zcolv = Z(:); % Make Z a column vector
    Const = ones(size(Xcolv)); % Vector of ones for constant term

    % find the coeffcients of the best plane fit
    Coefficients = [Xcolv Ycolv Const]\Zcolv; % Find the coefficients
    XCoeff = Coefficients(1); % X coefficient
    YCoeff = Coefficients(2); % X coefficient
    CCoeff = Coefficients(3); % constant term

    % detrend the data
    Z_p = XCoeff * X + YCoeff * Y + CCoeff;
    Z_f = Z - Z_p;

end