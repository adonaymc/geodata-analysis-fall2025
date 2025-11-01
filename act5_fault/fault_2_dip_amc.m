% CERI 7104/8104 Data Analysis in Geophysics
%
%
% - analyze stratigraphic layer, displaced by normal fault
% - determine dip of normal fault from cross-section
%
% @author thgoebel, CERI U of Memphis
% example from Trauth, Earth Sciences Recipes for Matlab, 3rd edtion
%
% Modified by Adonay Martinez-Coto 

clear all
clf
close all
%=======================1==================================================
%                  files and params
%==========================================================================
data_dir = '/home/pequenojn/gdrive/AdoMarti/01_Projects/DataAnalysis-Fall2025/homeworks/act5_fault/data/';
%'H:\Mi unidad\AdoMarti\01_Projects\DataAnalysis-Fall2025\homeworks\act5_fault\data\';
file_in  = "normalfault.txt"; % input data file 

rotAng   = 25;
%cmin   = -30;cmax = 22;
step    = 1;% desired grid step for interpolation
cd  % display current working directory
%=======================2==================================================
%                  load data, inspect data
%==========================================================================
XYZ = load(strcat(data_dir, file_in)); % load data

% Inspect data ranges
fprintf( '%s %.2f %.2f \n','X range:',min(XYZ(:,1)), max(XYZ(:,1)))
fprintf( '%s %.2f %.2f \n','Y range:',min(XYZ(:,2)), max(XYZ(:,2)))
fprintf( '%s %.2f %.2f \n','Z range (Topography):',min(XYZ(:,3)), max(XYZ(:,3)))

% move z to center
XYZ(:,1) = XYZ(:,1)-mean(XYZ(:,1));
XYZ(:,2) = XYZ(:,2)-mean(XYZ(:,2));

xmin = floor(min(XYZ(:,1))); xmax = ceil(max(XYZ(:,1)));
ymin = floor(min(XYZ(:,2))); ymax = ceil(max(XYZ(:,2)));
cmin = floor(min(XYZ(:,3))); cmax = ceil(max(XYZ(:,3)));

%=======================3==================================================
%                   rotation .   and interpolation
%==========================================================================
a_x = xmin:step:xmax;
a_y = ymin:step:ymax;
[XX, YY] = meshgrid( a_x, a_y); 
% % interpolation: nearest, linear, natural, cubic, v4 (biharmonic spline)
ZZ = griddata( XYZ(:,1), XYZ(:,2), XYZ(:,3), XX, YY, 'v4');

%plot 3D surface
figure(1)
set( gcf, 'Position', [1 1 2400 500]);
subplot(121)
mesh( XX, YY, ZZ), view( 0, 60)
colormap( 'hot'); 
cbar = colorbar;
cbar.Label.String = 'Z (m)';
xlabel( 'Easting (m)')
ylabel( 'Topography (m)')
xlim( [xmin, xmax])
ylim( [ymin, ymax])
zlim( [cmin, cmax])


%--------------rotate 3D matrix------------
[XYZ_rot] = rot_Mz(rotAng, XYZ);
ZZ_rot    = griddata( XYZ_rot(:,1), XYZ_rot(:,2), XYZ_rot(:,3), XX, YY, 'v4');

%plot rotated 3D surface
subplot(122)
mesh( XX, YY, ZZ_rot), view( 0, 60) %view(0,90)
colormap( 'hot');
cbar = colorbar;
cbar.Label.String = 'Z (m)';
xlabel( 'Easting (m)')
ylabel( 'Topography (m)')
xlim( [xmin, xmax])
ylim( [ymin, ymax])
zlim( [cmin, cmax])

saveas(gcf, 'fault_3D_rotated.png', 'png')

%=======================4==================================================
%                 stack profiles
%==========================================================================
[row,col] = size(ZZ_rot);
a_stack = zeros( 1, col);
figure(3)
% hide legend entries for individual profiles
for i=1:1:row
    a_stack = a_stack + ZZ_rot(i,:);
    plot(a_x, ZZ_rot(i,:), 'HandleVisibility','off'), hold on
    pause( .1)
end
plot(a_x, a_stack./row, 'r','LineWidth',3, 'DisplayName', 'Mean Profile (stacked)');
plot( a_x, mean(ZZ_rot), 'k--', 'DisplayName', 'Mean Profile (mean())');

xlabel( 'Easting (m)')
ylabel( 'Topography (m)')
legend show


[x,y] = ginput(2);
sprintf( '%s%.1f%s', "fault dip: ", points2angle(x,y))

saveas(gcf, 'fault_profiles.png', 'png');

%=========================================================================
%                 functions
%=========================================================================

function [angle] = points2angle(x,y)
    %
    % - compute angle between two points in cartesian
    %    using arctan(y2-y1, x2-x1)
    angle = rad2deg(atan2(y(2)-y(1),x(2)-x(1)));
end

function [M_rot] = rot_Mz(alpha, M)
    %
    %   rotates matrix M about z-axis
    %   using angle alpha
    %
    alpha = deg2rad(alpha);
    %rotation matrix
    %if str_axis == 'x'
    %rotx = [1 0 0; 0 cos(alpha) -sin(alpha) ; 0 sin(alpha) cos(alpha)] 
    %roty = [cos(alpha) 0 sin(alpha) ; 0 1 0 ; -sin(alpha) 0  cos(alpha)];
    rotz  = [cos(alpha) -sin(alpha) 0 ; sin(alpha) cos(alpha) 0 ; 0 0 1];
    M_rot = (rotz*M')';
end


