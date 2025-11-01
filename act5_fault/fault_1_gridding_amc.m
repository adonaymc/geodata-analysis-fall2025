% CERI 8104 / CIVL 8126 Data Analysis in Geophysics
%
% - analyze stratigraphic layer, displaced by normal fault
% - determine dip of normal fault from cross-section
%
% @author thgoebel, CERI U of Memphis
% example from Trauth, Earth Sciences Recipes for Matlab, 3rd edtion
%
% Modified by Adonay Martinez-Coto 

%=======================1==================================================
%                  files and params
%==========================================================================
data_dir = 'H:\Mi unidad\AdoMarti\01_Projects\DataAnalysis-Fall2025\homeworks\act5_fault\data\';
file_in  = "normalfault.txt"; % input data file 

step = 1; % desired grid step for interpolation
cd  % display current working directory

%=======================2==================================================
%                  load data, inspect data
%==========================================================================
XYZ = load(strcat(data_dir, file_in)); % load data

% here I decided to change the range min and max on the data and floor and ceil them to get integer values
xmin = floor(min(XYZ(:,1))); xmax = ceil(max(XYZ(:,1)));
ymin = floor(min(XYZ(:,2))); ymax = ceil(max(XYZ(:,2)));
cmin = floor(min(XYZ(:,3))); cmax = ceil(max(XYZ(:,3)));

% Inspect data ranges
fprintf( '%s %.2f %.2f \n','X range:',min(XYZ(:,1)), max(XYZ(:,1)))
fprintf( '%s %.2f %.2f \n','Y range:',min(XYZ(:,2)), max(XYZ(:,2)))
fprintf( '%s %.2f %.2f \n','Z range (Topography):',min(XYZ(:,3)), max(XYZ(:,3)))


labels = num2str( XYZ(:,3), 2);%vector with Z coordinate as string

% visualize raw data
fig1 = figure(1);
set( fig1, 'Position', [1 1 2400 500], 'Name', 'Normal Fault');
clf
subplot(121)
plot( XYZ(:,1),  XYZ(:,2), 'o'), hold on
text( XYZ(:,1)+.5, XYZ(:,2), labels, "FontSize", 7), hold on
xlim( [xmin, xmax])
ylim( [ymin, ymax])
xlabel( 'Easting (m)')
ylabel( 'Northing (m)')

%=======================3==================================================
%                  interpolation, gridding, surface plots
%==========================================================================
a_x = xmin:step:xmax;a_y = ymin:step:ymax;
[XX, YY] = meshgrid( a_x, a_y);% desired 
% interpolation: nearest, linear, natural, cubic, v4 (biharmonic spline)
ZZ = griddata( XYZ(:,1), XYZ(:,2), XYZ(:,3), XX, YY, 'v4');
%plot interpolation surface
contour( XX, YY, ZZ, 10)

% specify contour intervals, plot together with 2D cmap
%figure(2)
%clf
subplot(122)
a_con = cmin:10:cmax;
pcolor(XX, YY, ZZ), shading flat, hold on;
cbar = colorbar;
cbar.Label.String = 'Z (m)';
[c,h] = contour(XX,YY,ZZ, a_con, 'k');
clabel( c, h)
xlim( [xmin, xmax])
ylim( [ymin, ymax])
xlabel( 'Easting (m)')
ylabel( 'Northing (m)')

saveas(gcf, 'fault_2D_contour.png', 'png')

%----------plot 3D surface-------------
fig2 = figure(2);
set( fig2, 'Position', [1 1 2400 500], 'Name', 'Normal Fault - 3D Surface');
subplot(121)
% plot of a 3D meshed surface
mesh( XX, YY, ZZ), view(-37.5, 30) %view(0,90)
colormap( 'hot'); 
cbar = colorbar;
cbar.Label.String = 'Z (m)';
xlabel( 'Easting (m)')
ylabel( 'Northing (m)')
zlabel( 'Z (m)')
xlim( [xmin, xmax])
ylim( [ymin, ymax])
zlim( [cmin, cmax])

subplot(122)
% plot 3D contour surface
surf( XX,YY,ZZ), colormap( 'hot'); 
cbar = colorbar;
cbar.Label.String = 'Z (m)';
xlabel( 'Easting (m)')
ylabel( 'Northing (m)')
zlabel( 'Z (m)')
xlim( [xmin, xmax])
ylim( [ymin, ymax])
zlim( [cmin, cmax])

saveas(gcf, 'fault_3D_surface.png', 'png')