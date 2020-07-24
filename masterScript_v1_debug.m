function masterScript_v1_debug(h,pow,pp,dis,v0,operationmode,icecondition,vl,vr,scale,Ship,res,nLd,ice,i,rdag)


%the script is originally coded by Lauri Kuuliala. Fang Li modified this
%script by fixing some errors which could lead to unstable behavior,
%introducting crushing pattern on midship area to diminish the abnormal
%lateral force, and adding options for dynamic benging and random cusp
%size.

%the loop over i is for different cases and loop over j is for the segments
%of each case

%note that two versions of the code are needed for this, be sure to make
%any changes to both versions, they are different only in the start of the
%code during initialization

% icecondition: 1 for level ice; 2 for channel
% scale: scale factor for model scale test
% opetationmode: 1 for self-propulsion; 2 for towing
% h: ice thickness; std: stand. dev.
% pow: power in kW; pp: propeller pitch in percentage, 90 for 90%
% dis: distance of the simulation, m
% v0: initial speed
% vl, vr: left and right ice sheet drifting speed in converging ice channel
load DB10082018.mat
load areacorrec-h30.mat
warning('off')
% Default values:
% [operationmode,icecondition,vl,vr,scale] = [1, 1, 0, 0, 1];
if nargin == 6
    sigc = 1.28e6;
    stdsigc = 0.38e6;
    sigf = 0.404e6;
    stdsigf = 59.17e3;
    E = 5e9;
    operationmode = 1;
    icecondition = 1;
    vl = 0;
    vr = 0;
    scale = 1;
elseif nargin == 11
    operationmode = 1;
    icecondition = 1;
    vl = 0;
    vr = 0;
    scale = 1;    
end

tic
%background, make ship struct and hydrodynamic derivatives
% load Agulhashull.mat%hull line and angles
% x_f = x_f/scale;
% y_f = y_f/scale;
% x_phi = x_phi/scale;

vinit = v0;
%Thrust leave as [] if not known
Thrust = [];
power = pow*ones(1,dis);
prop_pitch = ones(1,dis)*pp*1.14/100;

%generise ship information
% [ship] = hull_teho(x_f, y_f, x_phi, phiIn, vinit, scale);
ship.name = 'Agu';
% 
% % Ice field thickness in 1m resolution
% iceThick = h+std*randn(1,dis);
% % eliminate abnormally large values
% for ii = 1:dis
%     while iceThick(ii)<h/4 || iceThick(ii)>h*7/4
%         iceThick(ii) = h+std*randn(1);
%     end
% end
% ridgeThick = zeros(length(iceThick),1);%Ice ridge thickness

%length of simulated time between intermediate saving
T = 3.000;
%time step size
dt = 0.002;
%time duration
t_total = inf;

lin = 0;

%for defining filename
model = lin;
parameters = h*100;

%initiate the simulation
    %calls simulation code
% disp('Warmup starts')
% [res, Ship, ice, ~] = simulation_start([0 3], ship, ...
%    0.005, iceThick, ridgeThick, lin, power, Thrust, prop_pitch, ...
%    DB,icecondition,scale,sigc,stdsigc,sigf,stdsigf,E,areacorr);

% %this part has to do with memory management and saving intermediate
% %results
% %    rowInd = find(nLoads(:,1)==0,1);
% %    nLd = nLoads(1:rowInd,:);
% 
% X = res(end,2) + Ship.x.*cos(res(end,4)) - Ship.y.*sin(res(end,4));
% Y = res(end,3) + Ship.x.*sin(res(end,4)) + Ship.y.*cos(res(end,4));
% 
% xlvl = min(X) - 5;
% yends = [min(Y-15) max(Y+15)];
% [ ~, ~, ii] = polyxpoly(ice.X, ice.Y, [xlvl xlvl], yends);
% if size(ii,1)==2
%    ice.X = [xlvl ice.X(ii(1,1):ii(2,1)) xlvl xlvl+500 xlvl+500 xlvl];
%    ice.Y = [yends(1)-100 ice.Y(ii(1,1):ii(2,1)) yends(2)+100 yends(2)+100 yends(1)-100 yends(1)-100];
% end
% 
% filename = strcat(ship.name,'-lin',num2str(model),'-h',num2str(parameters));
% 
% %further time segments
% res(1:end-2,:) = [];
% res(1,1) = 0;
% res(2,1) = dt;
% res(1:2,8) = vinit;
% res(1:2,9:end) = 0;
% nLd = [];
% 
% i = 1;
% filename = strcat(ship.name,'-lin',num2str(model),'-h',num2str(parameters),'-warmup');
% save(filename, 'ice', 'Ship')
% disp('Warmup ends, simulation starts')
while i < 1000
       TT = [i*T-T i*T];
%        todisplay = ['t=',num2str(TT(1))];
%        disp(todisplay)
       %debug
       %TT = [132.3,132.5];
       [res2, Ship, ice, nLoads] = simulation_continue(ice, res(end,:), ...
           TT, Ship, dt, power,Thrust,prop_pitch, lin, DB,vl,vr,icecondition, operationmode,scale,areacorr,rdag);

       res = vertcat(res,res2(2:end,:));
%        rowInd = find(nLoads(6:end,1)==0,1);
%        nLoads = nLoads(1:rowInd+5,:);
%        nLd1 = nLd;
       nLd = vertcat(nLd, nLoads);
       
       X = res(end,2) + Ship.x.*cos(res(end,4)) - Ship.y.*sin(res(end,4));
       Y = res(end,3) + Ship.x.*sin(res(end,4)) + Ship.y.*cos(res(end,4));
       
       xlvl = min(X) - 15;
       yends = [min(Y-15) max(Y+15)];
       [ ~,  ~, ii] = polyxpoly(ice.X, ice.Y, [xlvl xlvl], yends);
%        if size(ii,1)==2
%            ice.X = [xlvl ice.X(ii(1,1):ii(2,1)) xlvl xlvl+500 xlvl+500 xlvl];
%            ice.Y = [yends(1)-500 ice.Y(ii(1,1):ii(2,1)) yends(2)+500 yends(2)+500 yends(1)-500 yends(1)-500];
%        end
       
       if res(end,2)+75>=dis+200 || TT(2)>t_total
            i = 10000;
       end
       
       if mod(i,1) == 0 || res(end,2)+70>=dis+200
          
%           eval(['delete ',strcat(filename,'.mat')]);
          filename = strcat(ship.name,'-lin',num2str(model),'-h',num2str(parameters));
          save(filename, 'res', 'ice', 'Ship', 'nLd')
          
       end
       i = i + 1;
end

eval(['delete ',strcat(filename,'.mat')]);
filename = strcat(ship.name,'-lin',num2str(model),'-h',num2str(parameters));
save(filename, 'res', 'ice', 'Ship', 'nLd')

toc
