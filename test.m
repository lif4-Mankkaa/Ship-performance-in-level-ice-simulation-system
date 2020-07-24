theta = [0,0,0];
vn = 2;
dvn_over_dt = 0;
phi = 45/180*pi;
dphi_over_dt = 1;
d2phi_over_dt2 = 0;
rhoi = 900;
rhow = 1025;
hi = 0.3;
L = 1.5;
Lprime = 1.2;
dt = 0.001;
mu_n = 0.05;
mu_i = 0.05;
y = [0,0,0];
x = [0,0,0];
hcons = hi*(1-rhoi/rhow)/tan(phi);
h1 = hi*(1-rhoi/rhow);
h2 = hi*(rhoi/rhow);

theta_w = 120;
k1 = 2*tand(theta_w/2);
k2 = 2*Lprime*tand(theta_w/2)/(L-Lprime));
S = Lprime*k1*Lprime/2 + (L-Lprime)*k2*(L-Lprime)/2;
Sc = k1*(1/2*L*Lprime^2-1/3*Lprime^3) + 1/3*k2*(L-Lprime)^3;
SrhoBsq = k1*(1/2*L^2*Lprime^2-2/3*L*Lprime^3+1/4*Lprime^4) + 1/4*k2*(L-Lprime)^4;
m = [hi*S*rhoi;hi*S*rhoi;SrhoBsq*hi*rhoi-hi*S*rhoi*(Sc/S)^2];
% added_mass = [0.2,0.2,0.2];
% m = m.*(1+added_mass);

F = [];
the = [];
y_all = [];
x_all = [];
i = 1;
[theta_initial,s,vx,vy,stage] = rotation_initialize(vn,dvn_over_dt,phi,dphi_over_dt,d2phi_over_dt2,rhoi,rhow,hi,L,dt,mu_n,mu_i,S,Sc,SrhoBsq,hcons);
theta = theta_initial;
s_all(i) = s;
the(1:3,i) = theta;
y_all(1:3,i) = 0;
x_all(1:3,i) = 0;
while theta(1)<phi
    [F(1:3,i),theta_next,x_next,y_next,s_next,stage] = rotation_v2(m,theta,...
        vn,dvn_over_dt,phi,dphi_over_dt,d2phi_over_dt2,rhoi,rhow,hi,L,dt,mu_n,mu_i,S,Sc,SrhoBsq,y,x,s_all(i),hcons,stage);
    theta = theta_next;
    x = x_next;
    y = y_next;
    i = i+1;
    the(1:3,i) = theta;
    y_all(1:3,i) = y_next;
    x_all(1:3,i) = x_next;
    s_all(i) = s_next;
end
mean(F(1,:))
E_kinetic = 0.5*m(3)*theta_initial(2)^2+0.5*m(1)*(vx^2+vy^2);
E_potential = 0.25*rhow*9.81*Sc*Sc/S*(1-cos(2*theta_initial(1)));
E_initial = (E_kinetic+E_potential)*(sin(phi)+mu_n*cos(phi))^2
E_rotation = mean(F(1,:))*(sin(phi)+mu_n*cos(phi))*s_all(end)
vx_end = the(2,end)*Sc/S*sin(the(1,end))-x_all(2,end);
vy_end = y_all(2,end)-the(2,end)*Sc/S*cos(the(1,end));
vn_end = vx_end*sin(phi)-vy_end*cos(phi);
E_slamming = m(1)*(0-vn_end)*vn*sin(phi);
E_slamming = E_slamming*(sin(phi))^2
E_all = E_initial+E_rotation+E_slamming
% t_initial = s/vn
% theta_initial(1)/t_initial