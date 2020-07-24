function [ ship ] = hull_teho( x_frame, halfWidth, x_phi, phiIn, vinit, scale, T, zCOG)
%Definition of the data struct "ship".
%   Defines the discretised hull polygon based on frame spacing and half
%   widths. Flare angles are also determined based on a waterline below the
%   floating level. Unit tangents and (inward pointing) unit normals for
%   each discretisation point are also determined.

x_frame = x_frame(:);
halfWidth = halfWidth(:);
x_phi = x_phi(:);
phiIn = phiIn(:);

%delta_x for discretisation is minimum frame spacing / divisor
divisor = 400;
%//////////////////////////////////////////////////////////////////////////
spacing = abs(x_frame(1)-x_frame(end))/divisor;

%linear interpolation of points between frames for flotation wl
x_flot = x_frame(1):spacing:x_frame(end);
phi = interp1(x_phi,phiIn,x_flot)';
dphi_over_dx = -( phi(3:end)-phi(1:end-2) ) / ( 2*spacing );
dphi_over_dx = [dphi_over_dx(1); dphi_over_dx; 0];
d2phi_over_dx2 = ( phi(3:end)-2*phi(2:end-1)+phi(1:end-2) ) / spacing^2;
d2phi_over_dx2 = [d2phi_over_dx2(1); d2phi_over_dx2; 0];
y_flot = interp1(x_frame, halfWidth, x_flot);
x_flot = vertcat(x_flot(1:end)', flipud(x_flot(1:end)'));%oli äsken 1:end flipatussa
y_flot = vertcat(y_flot(1:end)', -1 * flipud(y_flot(1:end)'));
phi = vertcat(phi(1:end),flipud(phi(1:end)));
dphi_over_dx = vertcat(dphi_over_dx(1:end),flipud(dphi_over_dx(1:end)));
d2phi_over_dx2 = vertcat(d2phi_over_dx2(1:end),flipud(d2phi_over_dx2(1:end)));
%tangent is of the form [dx_tan; dy_tan]. the tangents are always in the
%negative x-direction
tangent = zeros(2,length(x_flot));
k = zeros(1, length(x_flot));
tangent(:,1) = [0 0];%OLI ENNEN [0 1]
for i = 2:length(x_flot)-1
    p = polyfit(x_flot(i-1:i+1), y_flot(i-1:i+1), 2);
    k(i) = 2*p(1)*x_flot(i) + p(2);
    if k(i) ~= 0 && isfinite(k(i))
        tangent(:,i) = [-1/sqrt(1+k(i)^2); -k(i)/sqrt(1+k(i)^2)];
    elseif k(i) == 0
        tangent(:,i) = [-1; 0];
    elseif isinf(k(i))
        tangent(:,i) = [0; 0];%OLI ENNEN [0 ;1]
    end
end
tangent(:,length(x_flot)) = [0 1];

%normal is of the form [dx_nor; dy_nor]. the normals point inside the hull
%polygon
normal = zeros(2, length(x_flot));
normal(:,1) = [1 0];
for i = 2:length(x_flot)-1
    if k(i) ~= 0 && isfinite(k(i))
        if y_flot(i) > 0
            normal(:,i) = [k(i)/sqrt(1+k(i)^2); -1/sqrt(1+k(i)^2)];
        else
            normal(:,i) = [-k(i)/sqrt(1+k(i)^2); 1/sqrt(1+k(i)^2)];
        end
    elseif k(i) == 0
        if y_flot(i) < 0
            normal(:,i) = [0; 1];
        else
            normal(:,i) = [0; -1];
        end
    elseif isinf(k(i))
        if x_flot(i) > 0
            normal(:,i) = [-1; 0];
        else
            normal(:,i) = [1; 0];
        end
    end
end

%directional angles of the tangents
dirAng = atan(k);

m = 12500e3/scale^3;
Ixx = 1.77e7/scale^5;
% Iyy = 1.77e10/scale^5;
Izz = 1.33e10/scale^5;

ship.x = single(x_flot);
ship.y = single(y_flot);
ship.phi = phi;
ship.dphi_over_dx = dphi_over_dx;
ship.d2phi_over_dx2 = d2phi_over_dx2;
ship.vcog = [vinit 0 0 0]';
ship.rcog = [-100/scale 100 0 0]';
ship.vow = 9/scale^0.5;
ship.M = [m 0 0 0;0 m 0 0;0 0 Ixx 0;0 0 0 Izz];
ship.F0 = [0 0 0 0]';
%Stability
ship.GMroll = 1.35/scale;
% ship.GMpitch = 45/scale;

ship.L = max(ship.x)-min(ship.x);
ship.B = max(ship.y) -min(ship.y);
ship.T = T/scale;
ship.zarm = (ship.T-zCOG)/scale; % vertical distance of gravity center to water surface
ship.Cb = (m/1024)/(ship.L*ship.B*ship.T);
ship.alpha = 33;
ship.tanalpha = tand(33);
ship.sinalpha = sind(33);
ship.cosalpha = cosd(33);
ship.phib = 20;
ship.tanphi = tand(20);
ship.sinphi = sind(20);
ship.cosphi = cosd(20);
ship.psi = atand(tand(ship.phib)/sind(ship.alpha));
ship.tanpsi = tand(ship.psi);
ship.sinpsi = sind(ship.psi);
ship.cospsi = cosd(ship.psi);
%ship.w = 0.2;%TÄM?IHAN HATUSTA
ship.delta = 0;

%ship sections
ship.max = max(ship.x);
ship.min = min(ship.x);
ind_mids = find(ship.y == ship.B/2,1)-1;
ship.mids = ship.x(ind_mids);
ind_mide = find(ship.y(ind_mids+1:end) < ship.B/2,1);
ship.mide = ship.x(ind_mide+ind_mids);
ship.Lbow = max(ship.x)-ship.mide; 
ship.stem = max(ship.x);

ship.normalx = single(normal(1,:));
ship.normaly = single(normal(2,:));
ship.normaly(length(ship.normaly)/2) = 0;%oli +0.5
ship.normalx(length(ship.normaly)/2) = -1;
ship.normaly(length(ship.normaly)/2+1) = 0;%oli +0.5
ship.normalx(length(ship.normaly)/2+1) = -1;
ship.normaly(end) = 0;
ship.normalx(end) = 1;
ship.normaly(1) = 0;
ship.normalx(1) = 1;
% ship.normaly(length(ship.normaly)/2+1.5) = 0;%oli +1,5
ship.tangentx = single(tangent(1,:));
ship.tangenty = single(tangent(2,:));
ship.tangenty(length(ship.tangenty)/2) = 0;
ship.tangentx(length(ship.tangenty)/2) = -1;
ship.tangenty(length(ship.tangenty)/2+1) = 0;
ship.tangentx(length(ship.tangenty)/2+1) = -1;
ship.tangenty(end) = 0;
ship.tangentx(end) = -1;
ship.tangenty(1) = 0;
ship.tangentx(1) = -1;
ship.dirAng = single(dirAng);

%hydrodynamic derivatives VERY CRUDELY ESTIMATED, USE BETTER DATA IF
%AVAILABLE

rhow = 1025;

const = -pi*(ship.T/ship.L)^2;
Yprimevdot = const*(1+0.16*ship.Cb*ship.B/ship.T-5.1*(ship.B/ship.L)^2);
Yprimerdot = const*(0.67*ship.B/ship.L-0.00033*(ship.B/ship.T)^2);
Yprimev = const*(1+0.4*ship.Cb*ship.B/ship.T);
Yprimer = const*(-0.5+2.2*ship.B/ship.L-0.08*ship.B/ship.T);
Nprimevdot = const*(1.1*ship.B/ship.L - 0.0003341*ship.B/ship.T);
Nprimerdot = const*(1/12 + 0.017*ship.Cb*ship.B/ship.T - 0.33*ship.B/ship.L);
Nprimev = const*(0.5 + 2.4*ship.T/ship.L);
Nprimer = const*(0.25 + 0.039*ship.B/ship.T - 0.56*ship.B/ship.L);

% % m11, m44 obtained assuming equivalent ellipsoid
% e = sqrt(1- (0.5*ship.B/2+0.5*ship.T)^2/(0.5*ship.L)^2);
% A0 = 2*(1-e^2)/e^3*(0.5*log((1+e)/(1-e))-e);
% B0 = 1/e^2-(1-e^2)/2/e^3*log((1+e)/(1-e));
% k11 = A0/(2-A0);
ship.A = 0.5*rhow*[0 0 0 0; 0 Yprimevdot*ship.L^3 0 Yprimerdot*ship.L^4; ...
    0 0 -2*0.007*3.14*ship.B^4*ship.L 0;0 Nprimevdot*ship.L^4 0 Nprimerdot*ship.L^5];
ship.Bp = 0.5*rhow*[0 0 0 0; 0 Yprimev*ship.L^2 0 Yprimer*ship.L^3; 0 0 0 0;0 Nprimev*ship.L^3 0 Nprimer*ship.L^4];

end

