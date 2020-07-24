function [F,theta_all_next,x_all_next,y_all_next,s_next,stage] = rotation_v3(m,theta,vn,dvn_over_dt,phi,dphi_over_dt,d2phi_over_dt2,rhoi,rhow,hi,L,dt,mu_n,mu_i,S,Sc,SrhoBsq,y,x,s,hcons,stage,itai)
% This is the function to calculate the force to turn an ice piece
% F: normal force onto ship hull, [Fx,Fy,Fz]
% theta_next: rotation angle at next time step in radian
% Lr: distance from center of rotation to origin at next step
% y_next: elevation of the far end from the beginning position at next step
% theta: rotation angle at current time step in radian
% vsinalpla_now: relative normal speed in horizontal plane at current time step
% phi: flare angle at current time step, in radian
% phi_next: flare angle at next time step, in radian
% rhoi, rhow: density of ice and water
% hi: ice thickness
% L: Length of ice piece
% dt: time step
% mu_n: friction coefficient in non-horizontal direction
% mu_i: ice-ice friction coefficient
% S: area of top surface
% Sc: first moment of ice piece top (or bottom) surface regarding far end
% SrhoBsq: 2nd moment of ice piece top (or bottom) surface regarding far end
% y_next: elevation of the far end from the biginning position at next step
hcons = 0;

dtheta_over_dt = theta(2);
theta = theta(1);
dx_over_dt = x(2);
x = x(1);
dy_over_dt = y(2);
y = y(1);

g = 9.81;
% theta0 = atan(hi/L);
% h1 = (1-rhoi/rhow)*hi;
% h2 = rhoi/rhow*hi;
h2 = 0;
c = Sc/S;
% cr = sqrt(hi/2^2+c^2);
cr = c;
% cr_left = sqrt(hi/2^2+(L-c)^2);
cr_left = L-c;
% thetac = atan(hi/2/c);
% thetac_left = atan(hi/2/(L-c));
thetac = 0;
thetac_left = 0;
% Lr: Diagonal length 
% Lr = sqrt(hi^2 + L^2);
Lr = L;
thetar = theta-thetac;
thetar_left = theta-thetac_left;
% m = [hi*S*rhoi;hi*S*rhoi;SrhoBsq*hi*rhoi-hi*S*rhoi*c^2];

if stage == 1
    % Calculate dtheta/dt, A*dtheta/dt = B
%     A = -(L*cos(theta)+L*sin(theta)*tan(phi));
%     B = -vn*tan(phi) + ( hcons+L*(1-cos(theta))-s )/(cos(phi))^2*dphi_over_dt;
%     dtheta_over_dt = B/A;
    A = -(L*cos(phi-theta));
    B = -vn*sin(phi) + (L-s)*(cos(phi))*dphi_over_dt - L*cos(phi-theta)*dphi_over_dt;
    dtheta_over_dt = B/A;   
    % Calculate d2theta/dt2, A*d2theta/dt2 = B
%     A = -(L*cos(theta)+L*sin(theta)*tan(phi));
%     B = ( L*cos(theta)*dtheta_over_dt^2-dvn_over_dt )*tan(phi)+...
%         2*( L*sin(theta)*dtheta_over_dt-vn )/cos(phi)^2*dphi_over_dt...
%         +( hcons+L*(1-cos(theta))-s )/cos(phi)^2*( d2phi_over_dt2+2*dphi_over_dt^2*tan(phi) )...
%         -L*sin(theta)*dtheta_over_dt^2;
    A = -L*cos(phi-theta);
    B = -2*vn*cos(phi)*dphi_over_dt - L*sin(theta)*cos(phi)/sin(phi)*d2phi_over_dt2 - ...
         dvn_over_dt*sin(phi) + L*sin(phi-theta)*dtheta_over_dt^2 - ...
         2*L*sin(phi-theta)*dtheta_over_dt*dphi_over_dt;
    d2theta_over_dt2 = B/A;
    
    % Matrix form: coef*F+const = m*accel
    % Acceleration
    accel = nan(3,1);
    accel(1) = d2theta_over_dt2*cr*sin(thetar)+dtheta_over_dt^2*cr*cos(thetar);
    accel(2) = -d2theta_over_dt2*cr*cos(thetar)+dtheta_over_dt^2*cr*sin(thetar);
    accel(3) = d2theta_over_dt2;
    % Force coefficient matrix
    coef = nan(3,3);
    coef(1,:) = [sin(phi)+mu_n*cos(phi),-1,0];
    coef(2,:) = [-cos(phi)+mu_n*sin(phi),0,-1];
    coef(3,:) = [cr_left*(cos(phi-thetar_left)-mu_n*sin(phi-thetar_left)),cr*sin(thetar),-cr*cos(thetar)];
    % Constants
    const = nan(3,1);
    const(1) = -rhow*g*(rhoi/rhow*hi-y+itai)*S*sin(theta)-rhow*g*Sc*sin(theta)^2;
    const(2) = -rhoi*g*hi*S+rhow*g*(rhoi/rhow*hi-y+itai)*S*cos(theta)+rhow*g*Sc*sin(theta)*cos(theta);
    const(3) = -rhow*g*sin(theta)*(SrhoBsq-S*c^2);
    F = coef\(m'.*accel-const);
    % Update theta
    % theta = atan(A)-atan(B)
    s_next = s + vn*dt;
%     A = -((hcons+L)*sin(phi)-s_next*sin(phi)+h2*cos(phi)-y*cos(phi))/Lr;
%     B = (L*tan(phi))/(-L);
%     theta_next = asin(A) - atan(B);
    theta_next = phi-asin( (L-s_next)*sin(phi)/L );
    
    if F(3)<0 || F(2)<0 || F(3)>mu_i*F(2)
        stage = 2;
    end
    
    theta_all_next = [theta_next,nan,nan];
    x_all_next = [0,0,nan];
    y_all_next = [0,0,nan];
end

if stage == 2
%     y = h2+L*sin(theta)-hi*cos(theta)+(hcons+L*(1-cos(theta))-hi*sin(theta)-s)*tan(phi);
    y = L*sin(theta) + (L*(1-cos(theta))-s)*tan(phi);
    dy_over_dt = (L*cos(theta)+hi*sin(theta)+L*sin(theta)*tan(phi)-hi*cos(theta)*tan(phi))*dtheta_over_dt...
        -vn*tan(phi) + ( hcons+L*(1-cos(theta))-hi*sin(theta)-s )/(cos(phi))^2*dphi_over_dt;
    coef = zeros(4,4);
    R = zeros(4,1);
    
    coef(1,:) = [sin(phi)+mu_n*cos(phi), -1, -m(1)*cr*sin(thetar), 0];
    coef(2,:) = [-cos(phi)+mu_n*sin(phi), -mu_i, m(2)*cr*cos(thetar), -m(2)];
    coef(3,:) = [cr_left*( cos(phi-thetar_left)-mu_n*sin(phi-thetar_left) ), cr*( sin(thetar)-mu_i*cos(thetar) ), -m(3), 0];
%     coef(4,:) = [0, 0, -( L*cos(theta)+L*sin(theta)*tan(phi) ), 1];
%     coef(4,:) = [0, 0, -L*cos(phi-theta)*cos(phi), cos(phi)^2];   
    coef(4,:) = [0, 0, -L*(sin(theta)+cos(theta)*cot(phi)),cot(phi)];
    
    R(1) = [rhow*g*(rhoi/rhow*hi-y+itai)*S*sin(theta)+rhow*g*Sc*sin(theta)^2+dtheta_over_dt^2*m(1)*cr*cos(thetar)];
    R(2) = [rhoi*g*hi*S-rhow*g*(rhoi/rhow*hi-y+itai)*S*cos(theta)-rhow*g*Sc*sin(theta)*cos(theta)+dtheta_over_dt^2*m(2)*cr*sin(thetar)];
    R(3) = [rhow*g*sin(theta)*(SrhoBsq-S*c^2)];
%     R(4) = [( L*cos(theta)*dtheta_over_dt^2-dvn_over_dt )*tan(phi)...
%         + 2*( L*sin(theta)*dtheta_over_dt-vn )/(cos(phi)^2)*dphi_over_dt ...
%         + ( hcons+L*(1-cos(theta))-s )/(cos(phi)^2)*( d2phi_over_dt2+2*dphi_over_dt^2*tan(phi) ) ...
%         - L*sin(theta)*dtheta_over_dt^2];
%     R(4) = L*sin(phi-theta)*cos(phi)*dtheta_over_dt^2 - dvn_over_dt*sin(phi)*cos(phi) ...
%         + 2*( L*sin(theta)*dtheta_over_dt-vn )*dphi_over_dt ...
%         + ( y-L*sin(theta) )/tan(phi)*d2phi_over_dt2 + 2*dphi_over_dt^2*( y-L*sin(theta) );
    R(4) = L*cos(theta)*dtheta_over_dt^2 - dvn_over_dt - ...
        (L*sin(theta)*dtheta_over_dt^2*cot(phi)-2*(dy_over_dt-L*cos(theta)*dtheta_over_dt)/sin(phi)^2*dphi_over_dt) ...
        + (y-L*sin(theta))/sin(phi)^2*(d2phi_over_dt2+2*dphi_over_dt^2*cot(phi));
    
    X = coef \ R;
    F = X(1:2);
    F(3) = X(2)*mu_i;
    d2theta_over_dt2 = X(3);
    d2y_over_dt2 = X(4);
    
    dtheta_over_dt_next = d2theta_over_dt2*dt + dtheta_over_dt;
    theta_next = theta+dtheta_over_dt*dt+0.5*d2theta_over_dt2*dt^2;
    s_next = s + vn*dt;
    dy_over_dt_next = d2y_over_dt2*dt + dy_over_dt;
    y_next = y + dy_over_dt*dt + 0.5*d2y_over_dt2*dt^2;
    
    if F(2)<0 || X(4)<0
        stage = 3;
    end
    
    theta_all_next = [theta_next,dtheta_over_dt_next,d2theta_over_dt2];
    x_all_next = [0,0,0];
    y_all_next = [y_next,dy_over_dt_next,d2y_over_dt2];
end

if stage == 3
%     y = h2+L*sin(theta)-hi*cos(theta)+(hcons+L*(1-cos(theta))-hi*sin(theta)-s)*tan(phi);
%     dy_over_dt = (L*cos(theta)+hi*sin(theta)+L*sin(theta)*tan(phi)-hi*cos(theta)*tan(phi))*dtheta_over_dt...
%         -vn*tan(phi) + ( hcons+L*(1-cos(theta))-hi*sin(theta)-s )/(cos(phi))^2*dphi_over_dt;
    coef = zeros(4,4);
    R = zeros(4,1);    

    coef(1,:) = [sin(phi)+mu_n*cos(phi), -m(1)*cr*sin(thetar), m(1), 0];
    coef(2,:) = [-cos(phi)+mu_n*sin(phi), m(2)*cr*cos(thetar), 0, -m(2)];
    coef(3,:) = [cr_left*( cos(phi-thetar_left)-mu_n*sin(phi-thetar_left) ), -m(3), 0, 0];
%     coef(4,:) = [0, -( L*cos(theta)+L*sin(theta)*tan(phi) ), tan(phi), 1];
%     coef(4,:) = [0, -L*cos(phi-theta)*cos(phi), sin(phi)*cos(phi), cos(phi)^2];
    coef(4,:) = [0, -L*(sin(theta)+cos(theta)*cot(phi)),1,cot(phi)];
    
    R(1) = [rhow*g*(rhoi/rhow*hi-y+itai)*S*sin(theta)+rhow*g*Sc*sin(theta)^2+dtheta_over_dt^2*m(1)*cr*cos(thetar)];
    R(2) = [rhoi*g*hi*S-rhow*g*(rhoi/rhow*hi-y+itai)*S*cos(theta)-rhow*g*Sc*sin(theta)*cos(theta)+dtheta_over_dt^2*m(2)*cr*sin(thetar)];
    R(3) = [rhow*g*sin(theta)*(SrhoBsq-S*c^2)];
%     R(4) = [( L*cos(theta)*dtheta_over_dt^2-dvn_over_dt )*tan(phi)...
%         + 2*( L*sin(theta)*dtheta_over_dt-vn-dx_over_dt )/(cos(phi)^2)*dphi_over_dt ...
%         + ( hcons+L*(1-cos(theta))-s-x )/(cos(phi)^2)*( d2phi_over_dt2+2*dphi_over_dt^2*tan(phi) ) ...
%         - L*sin(theta)*dtheta_over_dt^2];
%     R(4) = L*sin(phi-theta)*cos(phi)*dtheta_over_dt^2 - dvn_over_dt*sin(phi)*cos(phi) ...
%         + 2*( L*sin(theta)*dtheta_over_dt-vn-dx_over_dt )*dphi_over_dt ...
%         + ( y-L*sin(theta) )/tan(phi)*d2phi_over_dt2 + 2*dphi_over_dt^2*( y-L*sin(theta) );
    R(4) = L*cos(theta)*dtheta_over_dt^2 - dvn_over_dt - ...
        (L*sin(theta)*dtheta_over_dt^2*cot(phi)-2*(dy_over_dt-L*cos(theta)*dtheta_over_dt)/sin(phi)^2*dphi_over_dt) ...
        + (y-L*sin(theta))/sin(phi)^2*(d2phi_over_dt2+2*dphi_over_dt^2*cot(phi));
    
    X = coef \ R;
    F = X(1);
    F(2:3) = 0;
    d2theta_over_dt2 = X(2);
    d2x_over_dt2 = X(3);
    d2y_over_dt2 = X(4);
    
    dtheta_over_dt_next = d2theta_over_dt2*dt + dtheta_over_dt;
    theta_next = theta+dtheta_over_dt*dt+0.5*d2theta_over_dt2*dt^2;
    s_next = s + vn*dt;
    
    dx_over_dt_next = d2x_over_dt2*dt+dx_over_dt;
    x_next = dx_over_dt*dt + x + 0.5*d2x_over_dt2*dt^2;
    
    dy_over_dt_next = d2y_over_dt2*dt+dy_over_dt;
    y_next = dy_over_dt*dt + y + 0.5*d2y_over_dt2*dt^2;
    
    theta_all_next = [theta_next,dtheta_over_dt_next,d2theta_over_dt2];
    x_all_next = [x_next,dx_over_dt_next,d2x_over_dt2];
    y_all_next = [y_next,dy_over_dt_next,d2y_over_dt2];
end
% if dtheta_over_dt_next(1)<0
%     disp('Turning back')
% end
if isnan(F(1))
    disp('F_rotate is nan') 
end
if abs(F(1))>1e8
    disp('abnormal F_rotate')
end