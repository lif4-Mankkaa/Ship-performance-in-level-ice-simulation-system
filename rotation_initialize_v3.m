function [theta_initial,s,vx,vy,stage] = rotation_initialize_v3(vn,dvn_over_dt,phi,dphi_over_dt,d2phi_over_dt2,rhoi,rhow,hi,L,dt,mu_n,mu_i,S,Sc,SrhoBsq,hcons,itai)
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

% LET dx_tip/dt = dy_tip/dt = 0
theta = 0;
hcons = 0;

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
m = [hi*S*rhoi;hi*S*rhoi;SrhoBsq*hi*rhoi-hi*S*rhoi*c^2];
% Try stage 1
F = -1;
while F(1)<0
    theta = theta + 0.005;
    thetar = theta-thetac;
    thetar_left = theta-thetac_left;
%     s  = -(-L*sin(theta))/tan(phi)+L*(1-cos(theta));
    s = L - L*sin(phi-theta)/sin(phi);
%     A = -(L*cos(theta)+L*sin(theta)*tan(phi));
%     B = -vn*tan(phi) + ( L*(1-cos(theta))-s )/(cos(phi))^2*dphi_over_dt;
%     dtheta_over_dt = B/A;
    A = -(L*cos(phi-theta));
    B = -vn*sin(phi) + (L-s)*(cos(phi))*dphi_over_dt - L*cos(phi-theta)*dphi_over_dt;
    dtheta_over_dt = B/A;   
%     A = -(L*cos(theta)+L*sin(theta)*tan(phi));
%     B = ( L*cos(theta)*dtheta_over_dt^2-dvn_over_dt )*tan(phi)+...
%         2*( L*sin(theta)*dtheta_over_dt-vn )/cos(phi)^2*dphi_over_dt...
%         +( L*(1-cos(theta))-s )/cos(phi)^2*( d2phi_over_dt2+2*dphi_over_dt^2*tan(phi) )...
%         -L*sin(theta)*dtheta_over_dt^2;
%     d2theta_over_dt2 = B/A;
    A = -L*cos(phi-theta);
    B = -2*vn*cos(phi)*dphi_over_dt - L*sin(theta)*cos(phi)/sin(phi)*d2phi_over_dt2 - ...
         dvn_over_dt*sin(phi) + L*sin(phi-theta)*dtheta_over_dt^2 - ...
         2*L*sin(phi-theta)*dtheta_over_dt*dphi_over_dt;
    d2theta_over_dt2 = B/A;
    
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
    const(1) = -rhow*g*(rhoi/rhow*hi+itai)*S*sin(theta)-rhow*g*Sc*sin(theta)^2;
    const(2) = -rhoi*g*hi*S+rhow*g*(rhoi/rhow*hi+itai)*S*cos(theta)+rhow*g*Sc*sin(theta)*cos(theta);
    const(3) = -rhow*g*sin(theta)*(SrhoBsq-S*c^2);    
    
    F = coef\(m.*accel-const);
    
    if F(3)<0 || F(2)<0 || F(3)>mu_i*F(2) || theta>=phi
        stage = 2;
    else
        theta_initial = [theta, dtheta_over_dt, d2theta_over_dt2];
        stage = 1;
    end
end

if stage == 2
    theta = 0;
    F = -1;
    dy_over_dt = 0;
    while F(1)<0
        theta = theta + 0.001;
        thetar = theta-thetac;
        thetar_left = theta-thetac_left;
        s = L - L*sin(phi-theta)/sin(phi);
        y = 0;
%         A = -(L*cos(theta)+L*sin(theta)*tan(phi));
%         B = -vn*tan(phi) + ( hcons+L*(1-cos(theta))-hi*sin(theta)-s )/(cos(phi))^2*dphi_over_dt;
%         dtheta_over_dt = B/A;
        A = -(L*cos(phi-theta));
        B = -vn*sin(phi) + (L-s)*(cos(phi))*dphi_over_dt - L*cos(phi-theta)*dphi_over_dt;
        dtheta_over_dt = B/A; 
    %     dy_over_dt = (L*cos(theta)+hi*sin(theta)+L*sin(theta)*tan(phi)-hi*cos(theta)*tan(phi))*dtheta_over_dt...
    %         -vn*tan(phi) + ( hcons+L*(1-cos(theta))-hi*sin(theta)-s )/(cos(phi))^2*dphi_over_dt;
        coef = zeros(4,4);
        R = zeros(4,1);

        coef(1,:) = [sin(phi)+mu_n*cos(phi), -1, -m(1)*cr*sin(thetar), 0];
        coef(2,:) = [-cos(phi)+mu_n*sin(phi), -mu_i, m(2)*cr*cos(thetar), -m(2)];
        coef(3,:) = [cr_left*( cos(phi-thetar_left)-mu_n*sin(phi-thetar_left) ), cr*( sin(thetar)-mu_i*cos(thetar) ), -m(3), 0];
    %     coef(4,:) = [0, 0, -( L*cos(theta)+L*sin(theta)*tan(phi) ), 1];
%         coef(4,:) = [0, 0, -L*cos(phi-theta)*cos(phi), cos(phi)^2];
        coef(4,:) = [0, 0, -L*(sin(theta)+cos(theta)*cot(phi)),cot(phi)];

        R(1) = [rhow*g*(rhoi/rhow*hi-y+itai)*S*sin(theta)+rhow*g*Sc*sin(theta)^2+dtheta_over_dt^2*m(1)*cr*cos(thetar)];
        R(2) = [rhoi*g*hi*S-rhow*g*(rhoi/rhow*hi-y+itai)*S*cos(theta)-rhow*g*Sc*sin(theta)*cos(theta)+dtheta_over_dt^2*m(2)*cr*sin(thetar)];
        R(3) = [rhow*g*sin(theta)*(SrhoBsq-S*c^2)];
    %     R(4) = [( L*cos(theta)*dtheta_over_dt^2-dvn_over_dt )*tan(phi)...
    %         + 2*( L*sin(theta)*dtheta_over_dt-vn )/(cos(phi)^2)*dphi_over_dt ...
    %         + ( hcons+L*(1-cos(theta))-s )/(cos(phi)^2)*( d2phi_over_dt2+2*dphi_over_dt^2*tan(phi) ) ...
    %         - L*sin(theta)*dtheta_over_dt^2];
%         R(4) = L*sin(phi-theta)*cos(phi)*dtheta_over_dt^2 - dvn_over_dt*sin(phi)*cos(phi) ...
%             + 2*( L*sin(theta)*dtheta_over_dt-vn )*dphi_over_dt ...
%             + ( y-L*sin(theta) )/tan(phi)*d2phi_over_dt2 + 2*dphi_over_dt^2*( y-L*sin(theta) );
        R(4) = L*cos(theta)*dtheta_over_dt^2 - dvn_over_dt - ...
            (L*sin(theta)*dtheta_over_dt^2*cot(phi)-2*(dy_over_dt-L*cos(theta)*dtheta_over_dt)/sin(phi)^2*dphi_over_dt) ...
            + (y-L*sin(theta))/sin(phi)^2*(d2phi_over_dt2+2*dphi_over_dt^2*cot(phi));
    
        X = coef \ R;
        F = X(1:2);
        F(3) = X(2)*mu_i;
        d2theta_over_dt2 = X(3);
        d2y_over_dt2 = X(4);

        dtheta_over_dt_next = d2theta_over_dt2*dt + dtheta_over_dt;
        theta_next = theta+dtheta_over_dt*dt;
%         s_next = s + vn*dt;

        if F(2)<0 || theta>=phi || X(4)<0
            stage = 3;
        else
            theta_initial = [theta, dtheta_over_dt, d2theta_over_dt2];
            stage = 2;            
        end
    end
end

if stage == 3
    theta = 0;
    F = -1;
    while F(1)<0
        theta = theta + 0.005;
        thetar = theta-thetac;
        thetar_left = theta-thetac_left;
        s = L - L*sin(phi-theta)/sin(phi);
        y = 0;
        x = 0;
        A = -(L*cos(phi-theta));
        B = -vn*sin(phi) + (L-s)*(cos(phi))*dphi_over_dt - L*cos(phi-theta)*dphi_over_dt;
        dtheta_over_dt = B/A; 
       
        dx_over_dt = 0;
    %     dy_over_dt = (L*cos(theta)+hi*sin(theta)+L*sin(theta)*tan(phi)-hi*cos(theta)*tan(phi))*dtheta_over_dt...
    %         -vn*tan(phi) + ( hcons+L*(1-cos(theta))-hi*sin(theta)-s )/(cos(phi))^2*dphi_over_dt;
        coef = zeros(4,4);
        R = zeros(4,1);    

        coef(1,:) = [sin(phi)+mu_n*cos(phi), -m(1)*cr*sin(thetar), m(1), 0];
        coef(2,:) = [-cos(phi)+mu_n*sin(phi), m(2)*cr*cos(thetar), 0, -m(2)];
        coef(3,:) = [cr_left*( cos(phi-thetar_left)-mu_n*sin(phi-thetar_left) ), -m(3), 0, 0];
    %     coef(4,:) = [0, -( L*cos(theta)+L*sin(theta)*tan(phi) ), tan(phi), 1];
        coef(4,:) = [0, -L*cos(phi-theta)*cos(phi), sin(phi)*cos(phi), cos(phi)^2];
        coef(4,:) = [0, -L*(sin(theta)+cos(theta)*cot(phi)),1,cot(phi)];
        
        R(1) = [rhow*g*(rhoi/rhow*hi-y+itai)*S*sin(theta)+rhow*g*Sc*sin(theta)^2+dtheta_over_dt^2*m(1)*cr*cos(thetar)];
        R(2) = [rhoi*g*hi*S-rhow*g*(rhoi/rhow*hi-y+itai)*S*cos(theta)-rhow*g*Sc*sin(theta)*cos(theta)+dtheta_over_dt^2*m(2)*cr*sin(thetar)];
        R(3) = [rhow*g*sin(theta)*(SrhoBsq-S*c^2)];
    %     R(4) = [( L*cos(theta)*dtheta_over_dt^2-dvn_over_dt )*tan(phi)...
    %         + 2*( L*sin(theta)*dtheta_over_dt-vn-dx_over_dt )/(cos(phi)^2)*dphi_over_dt ...
    %         + ( hcons+L*(1-cos(theta))-s-x )/(cos(phi)^2)*( d2phi_over_dt2+2*dphi_over_dt^2*tan(phi) ) ...
    %         - L*sin(theta)*dtheta_over_dt^2];
        R(4) = L*sin(phi-theta)*cos(phi)*dtheta_over_dt^2 - dvn_over_dt*sin(phi)*cos(phi) ...
            + 2*( L*sin(theta)*dtheta_over_dt-vn-dx_over_dt )*dphi_over_dt ...
            + ( y-L*sin(theta) )/tan(phi)*d2phi_over_dt2 + 2*dphi_over_dt^2*( y-L*sin(theta) );   
        R(4) = L*cos(theta)*dtheta_over_dt^2 - dvn_over_dt - ...
            (L*sin(theta)*dtheta_over_dt^2*cot(phi)-2*(dy_over_dt-L*cos(theta)*dtheta_over_dt)/sin(phi)^2*dphi_over_dt) ...
            + (y-L*sin(theta))/sin(phi)^2*(d2phi_over_dt2+2*dphi_over_dt^2*cot(phi));
        
        X = coef \ R;
        F = X(1);
        F(2:3) = 0;
        d2theta_over_dt2 = X(2);
%         d2x_over_dt2 = X(3);
%         d2y_over_dt2 = X(4);
        if X(1)>0 && (X(3)<0)
            disp('x less than 0')
        end

        theta_initial = [theta, dtheta_over_dt, d2theta_over_dt2];
    end
end

vx = dtheta_over_dt*cr*sin(theta-thetac);
vy = -dtheta_over_dt*cr*cos(theta-thetac);
