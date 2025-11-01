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