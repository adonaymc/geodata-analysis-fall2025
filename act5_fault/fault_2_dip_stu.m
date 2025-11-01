% CERI 7104/8104 Data Analysis in Geophysics
%
%
% - analyze stratigraphic layer, displaced by normal fault
% - determine dip of normal fault from cross-section
%
% @author thgoebel, CERI U of Memphis
% example from Trauth, Earth Sciences Recipes for Matlab, 3rd edtion
clear all
clf
%=======================1==================================================
%                  files and params
%==========================================================================
data_dir = "???";
file_in  = "???";

rotAng   = 23;
%cmin   = -30;cmax = 22;
step    = ???;% desired grid step for interpolation
cd  % display current working directory
%=======================2==================================================
%                  load data, inspect data
%==========================================================================
XYZ = ???;
sprintf( '%s %.1f %.1f','range of topography:',min(XYZ(:,3)), max(XYZ(:,3)))

% move z to center
XYZ(:,1) = XYZ(:,1)-mean(XYZ(:,1));
XYZ(:,2) = XYZ(:,2)-mean(XYZ(:,2));

xmin = round(min(XYZ(:,1)))-5; xmax = round(max(XYZ(:,1)))+5;
ymin = round(min(XYZ(:,2)))-5; ymax = round(max(XYZ(:,2)))+5;

%=======================3==================================================
%                   rotation .   and interpolation
%==========================================================================
a_x = ???;
a_y = ???;
[XX, YY] = ???( a_x, a_y); 
% % interpolation: nearest, linear, natural, cubic, v4 (biharmonic spline)
ZZ = ???( XYZ(:,1), XYZ(:,2), XYZ(:,3), XX, YY, '???');
%plot 3D surface
figure(1)
???( ???, ???, ???), view( 0, 70)

%--------------rotate 3D matrix------------
[XYZ_rot] = rot_Mz(rotAng, XYZ);
ZZ_rot    = ???( XYZ_rot(:,1), XYZ_rot(:,2), XYZ_rot(:,3), XX, YY, 'v4');

figure(2)
mesh( XX, YY, ZZ_rot), view( 0, 70) %view(0,90)
%=======================5==================================================
%                 stack profiles
%==========================================================================
[row,col] = size(ZZ_rot)
a_stack = zeros( 1, col);
figure(3)
for i=1:1:row
    a_stack = a_stack + ZZ_rot(i,:);
    plot(a_x, ZZ_rot(i,:)), hold on
    pause( .1)
end
plot(a_x, a_stack./row, 'r','LineWidth',3);
plot( a_x, mean(ZZ_rot), 'k--')

xlabel( 'Easting (m)')
ylabel( 'Topography (m)')

[x,y] = ginput(2);
sprintf( '%s%.1f%s', "fault dip: ", points2angle(x,y))








