%
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
data_dir = ???;
file_in    = ???;

% if = 1 show pcolor of detrend surface
plot_surf = 1;
zMin = -40; % only for surface plots
zMax =  40;

%=======================2==================================================
%                   load data, gridding, detrending etc.
%==========================================================================
XYZ = load(strcat(???, ???));
%A% remove nans
XYZ =  XYZ(~isnan(XYZ(:, ???)),:);
X     =  XYZ(:, ???);
Y     =  XYZ(:, ???);
Z     =  XYZ(:, ???);

sprintf( '%s %.1f %.1f','range of topography:',min(Z), max(Z))
% spatial increments
fprintf( '%s%.2f', 'grid spacing: ', X(2) - X(1), ' m')
%B% interpolation to fill gaps from removed outliers
[XX, YY] = meshgrid( unique( X), unique( Y));
ZZ         = ???( ???, ???, ???, ???, ???, 'linear');
%D% detrend
ZZ_detr = detrend_2d( ???);

%C%visualize raw data
if plot_surf == 1
    figure(1)
    subplot(221)
    pcolor( ???, ???, ???); 
    title( 'raw data')
    shading interp;
    cbar = colorbar('eastoutside');
    cbar.Label.String = 'Roughness';
    
    subplot(222)
    title( 'detrended')
    pcolor( ???, ???, ???); 
    colormap hot
    shading flat;
    cbar = colorbar('eastoutside');
    caxis([zMin, zMax]);
    cbar.Label.String = 'Roughness';
    pause( 5)
end

% %=======================3==================================================
% %                                 stack PSDs of roughness profiles
% %==========================================================================
% [Ny,Nx] = size( ZZ);
% dx = XX(1,2)-XX(1,1);
% dy = YY(2,1)-YY(1,1);
% 
% % initialize matrizes for PSD profiles
% m_px=zeros( Ny, Nx);
% m_py=zeros( Nx, Ny);
% %-------stack all x -profiles
% for i=1:Ny-1
%     a_fftz     = ???;
%     m_py(i,:) = ???.^2*dy;
% end
% 
% for i=1:Nx-1
%     a_fftz      = ???;
%     m_px(i,:) = ???;
% end
% 
% % get the mean of all profiles
% a_Px=mean( ???);
% a_Py=mean( ???);
% 
% %=======================4==================================================
% %                              plot PSD and fit with ginput
% %==========================================================================
% a_freqY = ???;
% a_freqX = ???;
% figure(1)
% clf()
% %plot ave. power spectrum over wavelength
% loglog( ???, ???, 'r', 'LineWidth', 3);
% hold on
% loglog( ???, ???, 'b', 'LineWidth', 3);
% xlim( gca, [10, 5*1e3]);
% legend('Slip Perp.', 'Slip Parall.')
% ylabel('PSD (m^3)')
% xlabel('Wavelength (m)')
% 
% disp( 'Select Roughness Exp. Fitting Range')
% [x,y] = ginput(2);
% 
% %=======================4==================================================
% %                  fit roughness exponent
% %==========================================================================
% %-log-transform and lsq fit
% [par, R] = polyfit(  ???, ???, ???);
% gamma = par(1);
% Hurst    = (gamma(1)-1)*.5;
% sprintf(  "roughness exponent: %.1f, Hurst=%.2f", gamma, Hurst)







