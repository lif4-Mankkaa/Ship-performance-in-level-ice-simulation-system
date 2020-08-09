% Use this script to define ship and ice parameters and simulation scenario

%% Ship parameters
load Agulhashull22112019.mat % hull line and angles

T = 7.15; % Draught

alpha = 33; % Waterline angle at stem, deg
phib = 20; % Stem angle

zCOG = 8.88; % Location of CoG

m = 12500e3;
Ixx = 1.77e7;
Izz = 1.33e10;

[ship] = hull_teho(x_f, y_f, x_phi, phiIn, 1, T, zCOG, m, Ixx, Izz, alpha, phib);
ship.name = 'Agu'; % Name of the ship

% Propulsion
ship.dprop = 4.5; % propeller diameter, m
ship.fullpower = 9000; % Full power, kW
ship.fullpitch = 1.14; % Max propeller pitch to diameter ratio (p/D)
ship.vow = 9;
ship.GMroll = 1.35;

% A reference point of ship power, propeller p/D and speed
ship.refpower = 8820; % Reference power
ship.refpitch = 1.098; % reference p/D for bollard pull calculation
ship.refvow = 9; % Reference open water speed

% ship rudder
b = 4.5; % Rudder height
ship.Ar = 10.24; % Rudder area, m2
ship.Rasp = b^2/ship.Ar; % Aspect ratio
ship.chord = 2.5;% Chord length, m
ship.rudpos = [-57.5,4-zCOG]; % Rudder coordinate
ship.Kr = 0.5+0.5/(1+0.15/(2.4/ship.dprop));

% Round ship coordinate
ship.x = round(ship.x,4);
ship.y = round(ship.y,4);

%% Ice parameters
h = 0.3; % mean ice thickness, m
h_std = 0; % std of ice thickness, m

sigc = 1.28e6; % Mean compressive strength, Pa
sigc_std = 0.38e6; % Std of compressive strength, Pa
sigf = 4.04e5; % Mean flexural strength, Pa
sigf_std = 5.92e4; % Std of flexural strength, Pa
E = 5e9; % Young's modulus, Pa

%% Simulation scenario
rdag = 30; % Rudder angle, deg
pow = 4242; % Engine power, kW
pp = 57; % propeller pitch, % of max

dis = 2000; % Distance to simulate in ship initial heading direction, m
vinit = 4.3; % Initial speed, m/s
ship.vcog = [vinit 0 0 0]';

operationmode = 1; % 1 for self-propulsion, 2 for towing
icecondition = 1; % 1 for level ice, 2 for channel

% Compressive ice channel, not functioning in the current version
v_chan = [0,0]; 
vl = v_chan(1);
vr = v_chan(2);

Thrust = []; % Leave as [] if not known
power = pow*ones(1,dis);
prop_pitch = ones(1,dis)*pp*ship.fullpitch/100;

%% Ice field initialization

% Ice field thickness in 1m resolution
iceThick = h+h_std*randn(1,dis);
% eliminate abnormally large values
for i = 1:dis
    while iceThick(i)<0.75*h || iceThick(i)>1.25*h
        iceThick(i) = h+h_std*randn(1);
    end
end
ridgeThick = zeros(length(iceThick),1);%Ice ridge thickness

% Initiate ice field
m = 2000; % Dimension of ice sheet in y-direction
ship.rcog = [-100 m-200 0 0]';
n = 4000; % Dimension of ice sheet in x-direction
n = (length(iceThick)==1)*n+(length(iceThick)>1)*length(iceThick);%iceThick is scalar when mean value is used, vector when coordinates are used

d = 15; % Only useful in compressive channel
if icecondition == 1
    [iceField] = icefield(n,m, 0, iceThick,zeros(size(iceThick)),1,sigc,sigc_std,sigf,sigf_std,E); 
elseif icecondition == 2
    [iceField] = channelfield(n,m, d, iceThick, 1); 
end

%% Run
masterScript_v1(ship,iceField,h,dis,operationmode,icecondition,v_chan,1,power,prop_pitch,Thrust,vinit,rdag)