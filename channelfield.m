function [ iceField ] = channelfield( n,m,d,iceThick,scale)

%n is the length of the X-direction grid vector and m is the lenght of the
%Y-direction grid vector. Grid vector spacing is one (meters).
%Ice field properties are generated from normal distributions with means
%and standard deviations that are at the present tunable from this code.
    lwarmup = 0;
[iceField.N,iceField.M] = meshgrid(1:n+lwarmup, 1:m);

%Normally distributed:
if size(iceThick) == [1,1]
    %normal distribution and mean for ice thickness
    meanh = iceThick;
    stdevh = 0;%.01;
    lleg = 10;%leg length: lleg
    iceField.h = single(random('Normal', meanh, stdevh, m/lleg, n/lleg));
    for i = 1:m/lleg
        lengt = length(iceField.h(:,1));
        iceField.h = [iceField.h(1:lleg*(i-1)+1,:);ones(lleg-1,1)*iceField.h(lleg*(i-1)+1,:);iceField.h(lleg*(i-1)+2:lengt,:)];
    end
    for i = 1:n/lleg
        lengt = length(iceField.h(1,:));
        iceField.h = [iceField.h(:,1:lleg*(i-1)+1),iceField.h(:,lleg*(i-1)+1)*ones(1,lleg-1),iceField.h(:,lleg*(i-1)+2:lengt)];
    end
else
%\\\\\\\\\\\\\\old version v0\\\\\\\\\\\\\\\\\\
%     meanh = iceThick(1);
%     stedevh = 0;
%     iceField.thickness = iceThick;
%     iceField.h = ones(m,n)*iceThick(1);
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

    n = length(iceThick);
    meanh = [iceThick(1)*ones(lwarmup,1);smooth(iceThick,100)];
    warmupThick = mean(iceThick(1));%warmup length = 200
    stdevh = 0.01;
    iceField.h = [repmat(warmupThick*ones(1,lwarmup),m,1),repmat(iceThick,m,1)];
%     iceField.h = iceField.h + single(random('Normal', 0, stdevh, m, lwarmup+length(iceThick)));
end

% iceField.ridgethick = ridgethick;
%

%mean and stdev for sigmaf (ice flexural strength)
if scale == 1
    meansf = 0.404e6;
    stdevsf = 59.17e3;%.11e6;
    iceField.sigmaf = single(random('Normal', meansf, stdevsf, m, n+lwarmup)); 
    %mean and stdev for sigmac (ice crushing strength)
    meansc = 1.28e6;
    stdevsc = 0.38e6;%.5e6;
    iceField.sigmac = single(random('Normal', meansc, stdevsc, m, n+lwarmup));
    iceField.E = 5e9;
else
    meansf = 15e3;
    stdevsf = 3e3;%.11e6;
    iceField.sigmaf = single(random('Normal', meansf, stdevsf, m, n+lwarmup)); 
    %mean and stdev for sigmac (ice crushing strength)
    meansc = 12e3;
    stdevsc = 2e3;%.5e6;
    iceField.sigmac = single(random('Normal', meansc, stdevsc, m, n+lwarmup)); 
    iceField.E = 60e6;
end

%coefficient of friction, Poisson's ratio, density and Young's modulus of the ice field
iceField.mu = 0.1;
iceField.nu = 0.3;
% iceField.rho = 970;%[kg/m^3]
iceField.rhow = 1020;
iceField.rhoi = 900;%980;

%ice edge (now initially Y-axis straight line) 
paraXmean = 2.5;
paraXstd = 1;
lim1 = paraXmean-2*paraXstd;
lim2 = paraXmean+2*paraXstd;
paraYstd = 0.5;
rdnm = paraXstd*randn(1,10000)+paraXmean;
count = (rdnm<lim1) | (rdnm>lim2);
while sum(count) > 0
    nb = sum(count);
    rdnm((rdnm<lim1) | (rdnm>lim2)) = paraXstd*randn(1,nb)+paraXmean;
    count = (rdnm<lim1) | (rdnm>lim2);
end
    
Xedge1 = cumsum(rdnm);
Xedge1(Xedge1>n) = [];
Yedge1 = paraYstd*randn(1,length(Xedge1));

%repeat for the other edge
rdnm = paraXstd*randn(1,10000)+paraXmean;
count = (rdnm<lim1) | (rdnm>lim2);
while sum(count) > 0
    nb = sum(count);
    rdnm((rdnm<lim1) | (rdnm>lim2)) = paraXstd*randn(1,nb)+paraXmean;
    count = (rdnm<lim1) | (rdnm>lim2);
end
    
Xedge2 = cumsum(rdnm);
Xedge2(Xedge2>n) = [];
Yedge2 = paraYstd*randn(1,length(Xedge2));

% iceField.X1 = single([0*ones(1,(m-2*d)*5+1) Xedge1 n n 0]);
% iceField.Y1 = single([0:0.1:m/2-d m/2-d+Yedge1 m/2-d 0 0]);
% iceField.X2 = single([0 Xedge2 n n 0*ones(1,(m-2*d)*5+1)]);
% iceField.Y2 = single([m/2+d Yedge2+m/2+d m/2+d m m:-0.1:m/2+d]);
d = d/scale;
rightsec1X = 0*ones(1,(m-2*d)*5+1);
rightsec2X = [Xedge1,n]/scale;
leftsec1X = [n,Xedge2(end:-1:1)]/scale;
leftsec2X = 0*ones(1,(m-2*d)*5+1);
leftsec3X = (n+1)/scale;
rightsec3X = (n+1)/scale;
rightsec4X = 0;

rightsec1Y = 0:0.1:m/2-d;
rightsec2Y = [m/2-d+Yedge1/scale,m/2-d];
leftsec1Y = [m/2+d, Yedge2(end:-1:1)/scale+m/2+d];
leftsec2Y = m/2+d:0.1:m;
leftsec3Y = m;
rightsec3Y = 0;
rightsec4Y = 0;

iceField.X = single([rightsec1X rightsec2X leftsec1X ...
    leftsec2X leftsec3X rightsec3X rightsec4X]);
iceField.Y = single([rightsec1Y rightsec2Y leftsec1Y ...
    leftsec2Y leftsec3Y rightsec3Y rightsec4Y]);

iceField.dim = n/scale;

iceField.meanh = meanh;
iceField.meansf = meansf;
iceField.stdsf = stdevsf;
iceField.meansc = meansc;
iceField.stdsc = stdevsc;

% iceField.iceStart = single(iceStart);

end
