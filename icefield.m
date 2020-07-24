function [ iceField ] = icefield( n,m, iceStart, iceThick, ridgethick, scale,sigc,stdsigc,sigf,stdsigf,E)

%n is the length of the X-direction grid vector and m is the lenght of the
%Y-direction grid vector. Grid vector spacing is one (meters).
%Ice field properties are generated from normal distributions with means
%and standard deviations that are at the present tunable from this code.
lwarmup = 200;
[iceField.N,iceField.M] = meshgrid(1:n+lwarmup, 1:m);


%Normally distributed:
if size(iceThick) == [1,1]
    %normal distribution and mean for ice thickness
    meanh = iceThick;
    stdevh = 0.0935;%.01;
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
    iceField.h = [repmat(warmupThick*ones(1,lwarmup),m,1),repmat(iceThick,m,1)];
%     iceField.h = iceField.h + single(random('Normal', 0, stdevh, m, lwarmup+length(iceThick)));
end

iceField.ridgethick = ridgethick;
%

%mean and stdev for sigmaf (ice flexural strength)
if scale == 1
    meansf = sigf;%0.404e6;
    stdevsf = stdsigf;%59.17e3;%.11e6;
else
    meansf = sigf/scale;
    stdevsf = stdsigf/scale;
end
iceField.sigmaf = single(random('Normal', meansf, stdevsf, m, n+lwarmup)); 

%mean and stdev for sigmac (ice crushing strength)
if scale == 1
    meansc = sigc;%1.28e6;
    stdevsc = stdsigc;%0.38e6;%.5e6;
    iceField.E = E;
else
    meansc = sigc/scale;
    stdevsc = stdsigc/scale;
    iceField.E = E/scale;
end
iceField.sigmac = single(random('Normal', meansc, stdevsc, m, n+lwarmup));

%coefficient of friction, Poisson's ratio, density and Young's modulus of the ice field
iceField.mu = 0.1;
iceField.nu = 0.3;
% iceField.rho = 970;%[kg/m^3]
iceField.rhow = 1020;
iceField.rhoi = 900;%980;

%ice edge (now initially Y-axis straight line) 
iceField.X = single([iceStart*ones(1,m*10+1) max(n,length(iceThick))+lwarmup max(n,length(iceThick))+lwarmup iceStart]);
iceField.Y = single([(0:0.1:m) m 0 0]);
iceField.meanh = meanh;
iceField.meansf = meansf;
iceField.stdsf = stdevsf;
iceField.meansc = meansc;
iceField.stdsc = stdevsc;

iceField.iceStart = single(iceStart);
iceField.icepiece.num = 0;

iceField.sigmaf(1,1:5)
end
