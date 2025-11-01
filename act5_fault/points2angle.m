function [angle] = points2angle(x,y)
    %
    % - compute angle between two points in cartesian
    %    using arctan(y2-y1, x2-x1)
    angle = rad2deg(atan2(y(2)-y(1),x(2)-x(1)));
end