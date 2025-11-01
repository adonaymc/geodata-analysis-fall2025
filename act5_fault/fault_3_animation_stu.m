% CERI 7104/8104 Data Analysis in Geophysics
%
%
%           - rotate surface about z-axis, save figures
%               create animation
%
%       animation with ImageMagick:
% 1) convert   -delay 30   -loop 0   NormalFault_*.png   fault_rotate.gif
% 2) convert fault_rotate.gif -resize 300x300 fault_rotate_small.gif 
%   see: http://www.imagemagick.org/Usage/anim_basics/
%
% @author thgoebel, CERI U of Memphis
% 
clear all

%=======================1==================================================
%                  files and params
%==========================================================================
data_dir = "???";
file_in  = "???";

a_rotAng   = 0:20:360;
step       = 1;
%=======================2==================================================
%                  load data, inspect data
%==========================================================================
XYZ = ???;
sprintf( '%s %.1f %.1f','range of topography:',min(XYZ(:,3)), max(XYZ(:,3)))

% move z to center
XYZ(:,1) = ???(:,1)-???(XYZ(:,1));
XYZ(:,2) = ???(:,2)-???(XYZ(:,2));

xmin  = round(min(XYZ(:,1)))-5; 
xmax = round(max(XYZ(:,1)))+5;
ymin  = round(min(XYZ(:,2)))-5; 
ymax = round(max(XYZ(:,2)))+5;

%=======================3==================================================
%                    rotation .   and interpolation
%==========================================================================
a_x = xmin:step:xmax;
a_y = ymin:step:ymax;
[XX, YY] = ???( ???, ???); 

%--------------rotate 3D matrix------------
for i=1:length(a_rotAng)
    rotAng     = a_rotAng(i);
    file_out   = sprintf('animation/NormalFault_%03d.png', rotAng)
    [XYZ_rot]  = rot_Mz( ???, ???);

    ZZ_rot = ???( XYZ_rot(:,1), XYZ_rot(:,2), XYZ_rot(:,3), XX, YY, 'v4');

    fig = figure(1);
    clf
    % plot 3D contour surface of rotated, raster data
    ???( ???, ???, ???)
    saveas( fig, file_out)
end









