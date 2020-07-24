% Visualize

% theta = [0,0,0];
% vn = 1;
% dvn_over_dt = 0;
% phi = 45/180*pi;
% dphi_over_dt = 0;
% d2phi_over_dt2 = 0;
% rhoi = 900;
% rhow = 1025;
% hi = 0.3;
% L = 1.5;
% Lprime = 1.2;
% dt = 0.001;
% mu_n = 0.05;
% mu_i = 0.05;
% y = [0,0,0];
% x = [0,0,0];
% s = 0;
% stage = 1;
% h1 = hi*(1-rhoi/rhow);
% h2 = hi*(rhoi/rhow);
% hcons = h1/tan(phi);

% theta_w = 120;
% k1 = 2*tand(theta_w/2);
% k2 = 2*tand(Lprime*tand(theta_w/2))/(L-Lprime);
% S = Lprime*k1*Lprime/2 + (L-Lprime)*k2*(L-Lprime)/2;
% Sc = k1*(1/2*L*Lprime^2-1/3*Lprime^3) + 1/3*k2*(L-Lprime)^3;
% SrhoBsq = k1*(1/2*L^2*Lprime^2-2/3*L*Lprime^3+1/4*Lprime^4) + 1/4*k2*(L-Lprime)^4;
figure(1)
t = [1:length(F(1,:))]*0.001;
for i = 1:length(F(1,:))

    theta = the(1,i);
%     x_F = vn*t(i);
    x_F = s_all(i);
    y_F = 0;
    x_A = hcons+L*(1-cos(theta))-x_all(1,i);
    y_A = -h2+y_all(1,i)-L*sin(theta);
    x_B = x_A-hi*sin(theta);
    y_B = y_A+hi*cos(theta);
    x_Q = hcons+L-x_all(1,i);
    y_Q = -h2+y_all(1,i);
    x_Qup = x_Q-hi*sin(theta);
    y_Qup = y_Q+hi*cos(theta);
    
    fill([x_F-5,x_F-2,x_F,x_F+2,x_F-5,x_F-5],[-2*tan(phi),-2*tan(phi),0,2*tan(phi),2*tan(phi),-2*tan(phi)],[0.4 0.4 0.4])
    hold on
    fill([x_F-2,x_B,x_A,x_Q,L+3+hcons,L+3+hcons],[-2*tan(phi),y_B,y_A,0,0,-2*tan(phi)],'b')
    fill([x_A,x_B,x_Qup,x_Q,x_A],[y_A,y_B,y_Qup,y_Q,y_A],[0.9 0.9 0.9])
    line([-2,3],[0,0],'LineStyle',':')
    fill([L,L,L+3,L+3,L]+hcons,[h1,-h2,-h2,h1,h1],[0.9 0.9 0.9])
    hold off
    xlim([-2,3])
    ylim([-2,1])
    pause(0.001)
    
end