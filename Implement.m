% h = 0.23; % m
% h_std = 0.04; % m
% pow = 4500; % kW
% pp = 57; % percent of max


caseno = 1
if caseno == 2
    h = 0.22; % m
    h_std = 0; % m
    rdag = 30;
% pow = 4465; % kW
    pow = 4242; % kW
    pp = 57; % percent of max
elseif caseno == 1
    h = 0.29; % m
    h_std = 0; % m
    rdag = 28.5;
    pow = 4465; % kW
    pp = 70; % percent of max
end
dis = 2000; % m
v0 = 4.3; % m/s
sigc = 1.28e6; % Pa
sigc_std = 0.38e6; % Pa
sigf = 4.04e5; % Pa
sigf_std = 5.92e4; % Pa
E = 5e9; % Pa
operationmode = 1; % 1 for self-propulsion, 2 for towing
icecondition = 1; % 1 for level ice, 2 for channel
v_chan = [0,0];
vl = v_chan(1);
vr = v_chan(2);
scale = 1;


%% Ship and Ice field initialize
load Agulhashull22112019.mat%hull line and angles
x_f = x_f/scale;
y_f = y_f/scale;
x_phi = x_phi/scale;

vinit = v0;
%Thrust leave as [] if not known
Thrust = [];
power = pow*ones(1,dis);
prop_pitch = ones(1,dis)*pp*1.14/100;

%generise ship information
T = 7.15;
zCOG = 8.88;
[ship] = hull_teho(x_f, y_f, x_phi, phiIn, vinit, scale, T, zCOG);
ship.name = 'Agu';
%propeller
%propeller diameter, not scaled yet
ship.dprop = 4.5;
ship.fullpower = 8820;
ship.fullpitch = 1.14;
ship.refpitch = 1.098;%reference pitch for bollard pull calculation
ship.itaref = -0.330*(ship.refpitch)^2+0.958*ship.refpitch +0.018;
%ship rudder information, not scaled yet
b = 4.5;% Rudder height
ship.Ar = 10.24;
ship.Rasp = b^2/ship.Ar;
ship.chord = 2.5;%chord length
ship.rudpos = [57.5,zCOG-4]; % distance to CoG, distance to WL
ship.Kr = 0.5+0.5/(1+0.15/(2.4/ship.dprop));
%assumed 9MW propulsive power and 4.5m diameter propeller, not scaled yet
K_ref = -0.5848*ship.refpitch+3.349;
Ke_ref = (0.250/pi^2)^(1/3)*K_ref;
ship.bollard = 2*Ke_ref*(1/2*ship.fullpower*ship.dprop)^(2/3)*1e3;
ship.vow = 0.6595*(ship.fullpower*ship.itaref)^(0.3064);

ship.x = round(ship.x,4);
ship.y = round(ship.y,4);

% Ice field thickness in 1m resolution
iceThick = h+h_std*randn(1,dis);
% eliminate abnormally large values
for i = 1:dis
    while iceThick(i)<0.75*h || iceThick(i)>1.25*h
        iceThick(i) = h+h_std*randn(1);
    end
end
ridgeThick = zeros(length(iceThick),1);%Ice ridge thickness

m = 2000;
ship.rcog = [-100/scale m-200 0 0]';
n = 4000;
d = 15;
n = (length(iceThick)==1)*n+(length(iceThick)>1)*length(iceThick);%iceThick is scalar when mean value is used, vector when coordinates are used
if icecondition == 1
    [iceField] = icefield(n,m, 0, iceThick,zeros(size(iceThick)),scale,sigc,sigc_std,sigf,sigf_std,E); 
elseif icecondition == 2
    [iceField] = channelfield(n,m, d, iceThick, scale); 
end

%% Run
masterScript_v1(ship,iceField,h,dis,operationmode,icecondition,v_chan,scale,power,prop_pitch,Thrust,vinit,rdag)

% % % %% Debug
%  i = 160;
%  masterScript_v1_debug(h,pow,pp,dis,v0,operationmode,icecondition,vl,vr,scale,Ship,res,nLd,ice,i,rdag)