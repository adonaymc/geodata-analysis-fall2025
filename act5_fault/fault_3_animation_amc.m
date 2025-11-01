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
% example from Trauth, Earth Sciences Recipes for Matlab, 3rd edtion
%
% Modified by Adonay Martinez-Coto 



clear all
clf
close all

%=======================1==================================================
%                  files and params
%==========================================================================
data_dir = 'H:\Mi unidad\AdoMarti\01_Projects\DataAnalysis-Fall2025\homeworks\act5_fault\data\';
file_in  = "normalfault.txt"; % input data file 

a_rotAng   = 0:10:360;
step       = 1;

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
%                    rotation .   and interpolation
%==========================================================================
a_x = xmin:step:xmax;
a_y = ymin:step:ymax;
[XX, YY] = meshgrid( a_x, a_y);

fig = figure(1);
set( fig, 'Position', [1 1 700 600], 'Color', 'w');


% Here create gif and mp4 video
gif_filename = 'fault_rotate.gif';

video = VideoWriter('fault_rotate.mp4','MPEG-4');
video.FrameRate = 3; % Set the frame rate (frames per second)
open(video);

for i=1:length(a_rotAng)
    rotAng     = a_rotAng(i);
    file_out   = sprintf('animation/NormalFault_%03d.png', rotAng);
    [XYZ_rot]  = rot_Mz( rotAng, XYZ);

    ZZ_rot = griddata( XYZ_rot(:,1), XYZ_rot(:,2), XYZ_rot(:,3), XX, YY, 'v4');

    surfc(XX, YY, ZZ_rot);
    title(['Rotation Angle: ', num2str(rotAng), '°']);
    xlabel('Easting (m)');
    ylabel('Northing (m)');
    zlabel('Topography (m)');
    zlim([cmin cmax]);
    colormap('hot');
    colorbar;
    cbar = colorbar;
    cbar.Label.String = 'Z (m)';

    drawnow; % Force MATLAB to draw the plot now

    frame = getframe(fig); % Capture the current figure as a frame
    im = frame2im(frame); % Convert the frame to an rgb image
    % Convert to indexed image, gif are saved as an indexed image and a color map
    [imind, cm] = rgb2ind(im, 256);

    % Write to the GIF File
    if i == 1
        imwrite(imind, cm, gif_filename, 'gif',...
        'Loopcount', inf,... % create an infinite loop, it wont stop
        'DelayTime', 0.3); % DelayTime in seconds between frames.
    else
        imwrite(imind, cm, gif_filename, 'gif',...
        'WriteMode', 'append',... % append to the existing file
        'DelayTime', 0.5);
        % 
    end

    writeVideo(video, frame); % Write the frame to the video file

    pause(0.1);
end

close(video); % Close the video file

%==========================================================================
%                 functions
%==========================================================================

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






