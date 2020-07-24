function [ IN ] = inhull( xh, yh, xi, yi, xmax, xmin, ymax, ymin, xmid )
%detects ice points inside hull polygon
%   hull polygon is assumed to be convex with parallel midship and blunt
%   bow (two points with equal x at bow)
%   xh, yh are column vectors containing hull vertices
%   xi, yi are ice points in column vectors
%   xmax, ymax, xmin, ymin are minimum and maximum of hull polygon
%   xmid is 2*1 vector with start and end of paralle midship
%   in is vector with same size as xi and yi with 1 for points inside hull
%   and 0 for outside

%initialize IN
IN = false(size(xi));

%midship
IN(xi>=xmid(1)&xi<=xmid(2)&yi>ymin&yi<ymax) = true;
%bow y>0
IN(xi>xmid(2)&xi<xmax&yi>0&yi<interp1(xh(xh>=xmid(2)&yh>=0),yh(xh>=xmid(2)&yh>=0),xi)) = true;
%bow y<0
IN(xi>xmid(2)&xi<xmax&yi<0&yi>interp1(xh(xh>=xmid(2)&yh<=0),yh(xh>=xmid(2)&yh<=0),xi)) = true;
%stern y>0
IN(xi<xmid(1)&xi>xmin&yi>0&yi<interp1(xh(xh<=xmid(1)&yh>=0),yh(xh<=xmid(1)&yh>=0),xi)) = true;
%stern y<0
IN(xi<xmid(1)&xi>xmin&yi<0&yi>interp1(xh(xh<=xmid(1)&yh<=0),yh(xh<=xmid(1)&yh<=0),xi)) = true;
%y==0
IN(xi>xmin&xi<xmax&yi==0) = true;

end

%xh = ship.x;yh = ship.y;xi = iceField.x;yi = iceField.y;xmax =max(ship.x);
%ymax = max(ship.y),xmid = [-31.4918 26.5030];ymin=min(ship.y)
