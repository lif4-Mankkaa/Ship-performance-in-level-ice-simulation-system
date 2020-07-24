function [results, ship, iceField, normalLoads ] = ...
    simulation_continue(iceField,initial,tspan,ship,dt,power,Thrust,prop_pitch,lin,DB,vl,vr,icecondition,operationmode,scale,areacorr,rdag)
warning('off','all')
%Numerical simulation of ship operating in a level ice field
%definition of structs 'ship' and 'iceField'
%[ship, YDer, NDer] = hull_teho(x_flot, y_flot, x_bel, y_bel);
%size of the ice field in meters. m is the extent in the X- direction and n
%in Y- direction
% m = 30000;
% n = 1000;
% iceStrt = 170;
% [iceField, meanh] = icefield(m,n, iceStrt); 

if length(iceField.meanh) == 1
    meanh = iceField.meanh;
else
    meanh = iceField.meanh(floor(max(initial(2),1)));%160 is ice starting coordinate
end

ship.rcog = initial(2:5)';
ship.vcog = initial(6:9)';
ship.F0 = initial(14:17)';

%start and end time and time step
t = tspan(1);
tend = tspan(2);

%initial conditions
rk = ship.rcog;
vk = ship.vcog;
Fk = ship.F0;

%mass matrix, added masses taken into account in the slow motion
%derivatives
M = ship.M;
A = -ship.A;
BM = -ship.Bp;

ak = (M+A)^-1*ship.F0;
if operationmode == 2
    ak(1) = 0;
end

%output initialization. counter is the row of next set of results,
%nLnext the row of next normal loads
results = zeros(round((tend-t)/dt), 23);
results(1,:) = [t rk(1) rk(2) rk(3) rk(4) vk(1) vk(2) vk(3) vk(4) ak(1) ak(2) ak(3) ak(4) 0 0 0 0 0 0 meanh 0 0 0];
counter = 2;
normalLoads = zeros(5*round((tend-t))/dt, 17);
nLnext = 1;

rhow = iceField.rhow;
rhoi = iceField.rhoi;
mu = iceField.mu;

L = ship.L;
B = ship.B;
T = ship.T;

tanphi = ship.tanphi;
tanalpha = ship.tanalpha;
cosphi = ship.cosphi;
cospsi = ship.cospsi;
sinphi = ship.sinphi;

E = iceField.E;
nu = iceField.nu;
hi = meanh;

l_c = (E*hi.^3/(12*(1-nu.^2)*rhow*9.81)).^0.25;

%submersion force (this is constant)
Fbow_buo = (rhow - rhoi)*9.81*meanh*B*T*(B+T)/(B+2*T);
Fbow_fric = (rhow - rhoi)*9.81*meanh*B*mu*(0.7*L - T/tanphi - B/(4*tanalpha));
Fmid = (rhow - rhoi)*9.81*meanh*B*mu*T*cosphi*cospsi*sqrt(1/sinphi^2+1/tanalpha^2);
Rsub = Fbow_buo + Fbow_fric + Fmid;
% Rsub = (rhow - rhoi)*9.81*meanh*B*(T*(B+T)/(B+2*T) + mu*(0.7*L - T/tanphi - B/(4*tanalpha) + T*cosphi*cospsi*sqrt(1/sinphi^2+1/tanalpha^2)));

rudderangle = rdag;
while t < tend-0.0001 && rk(1) + ship.L < length(iceField.h(1,:))+ship.L*3
    psi0 = rk(4);
    TM = eye(4);
    TM(1:2,1:2) = [cos(psi0) -sin(psi0); sin(psi0) cos(psi0)]; 
    ri0 = rk + TM*(vk*dt + 0.5*ak*dt^2);
    vi0 = vk + ak*dt;
    rudderang = rudderangle - atand( (vi0(2)-vi0(4)*ship.rudpos(1))/vi0(1) );
    % ----------------------------------------------------------------------------------------------- %
    % Ice forces
    % Ridge resistance by Malmberg's method
    Rridge = 0;
    % Rridge = midshipResistance(T,B,iceField.ridgethick(tridge),69) + bowResistance(T,B,ship.psi,ship.alpha,mu,iceField.ridgethick(tridge));
    % Icebreaking
    [Fi0, iceField_t, nLoad] = iceForces(ri0, vi0, ship, iceField, [vl,vr], t, lin, DB, areacorr);
    % Ice rotate
    g = 9.81;
    Hbw = 3*(ship.B/ship.Lbow)*vi0(1)^2/(2*g);
    xc = 0.372*vi0(1)^2/g;
    k = g/vi0(1)^2;
    [F_rotate, iceField_t] = rotateForces(iceField_t, ship, ri0, vi0,rhoi,rhow,dt,mu,t,Hbw,xc,k);
    if ~isempty(Thrust)
        Ft = Thrust(min(round(t)+1,length(Thrust)));
    else
        power_t = power(min(round(t)+1,length(power)));
        pp_D = prop_pitch(min(round(t)+1,length(power)));
        Ft = propellerForce(vi0, ship, pp_D, power_t,85);
    end
    
    if operationmode == 2
        Ft(1) = 0;
    end
    
    Res = Rsub * (1 + 9.4 * norm(vi0(1:2))/sqrt(9.81*ship.L))+Rridge;
    % Assume acting point: 0.25L ahead of CoG; 0.5T below CoG
    v_action = [vi0(1),vi0(2)+vi0(4)*0.25*ship.L];
    Fxs = -Res * v_action(1)/sqrt((v_action(1))^2+(v_action(2))^2);
    Fxy = -Res * v_action(2)/sqrt((v_action(1))^2+(v_action(2))^2);
%     Fxz = (Fmid/mu+Fbow_buo/tand(ship.phib*2)) * (1 + 9.4 * norm(vi0(1:2))/sqrt(9.81*ship.L));
    Rs = [Fxs Fxy Fxy*0.5*ship.T Fxy*0.25*ship.L]';
    %rudder force
    Fr = 2*rudderforce(ship.Rasp,rudderang,ship.chord,vi0(1),ship.Ar,-ship.rudpos, 0.15, ship.Kr, Ft(1)/2, ship.dprop);
    F_nonlin = [M(1,1)*vi0(2)*vi0(4),-M(2,2)*vi0(1)*vi0(4),0,0]';
    Fboy = [0 0 -ship.GMroll*ri0(3)*ship.M(1,1)*9.8 0]';
    F0 = Fi0 + Ft + Rs + Fr + Fboy + F_rotate + F_nonlin;
    ai = (M+A)^-1*(F0-(vi0(1)*BM)*vi0);
    
    ri = ri0; vi = vi0; Fii = Fi0; Fi = F0;
    
    iceField = iceField_t;
    
    %         Ice field drift
    %            Find breakpoint 1 (separate channel edges), breakpoint 2 (separ. n,m and n,0)
    if icecondition == 2
        bk1 = find(iceField.X == min(iceField.X));
        bk1 = round(length(bk1/2));
        bk2 = find(iceField.X == iceField.dim);
        
        tempind = zeros(1,length(iceField.Y));
        tempind(bk1:bk2(1)) = 1;
        %             on the right hand side
        iceField.Y(tempind == 1) = iceField.Y(tempind == 1) + vl*dt;
        tempind = zeros(1,length(iceField.Y));
        tempind(bk2(2):end-bk1-2) = 1;
        %             on the left hand side
        iceField.Y(tempind == 1) = iceField.Y(tempind == 1) - vr*dt;
    end
    
    rk = ri;
    vk = vi;
    ak = ai;
%     if ai(2)<-0.0023
%         disp('a large')
%     end
    if operationmode == 2
        ak(1) = 0;
    end
    
    if ~isempty(nLoad)
        nLend = nLnext + length(nLoad(:,1)) -1;
        normalLoads(nLnext:nLend,1) = t*ones(length(nLoad(:,1)),1);
        normalLoads(nLnext:nLend,2:17) = nLoad;
        nLnext = nLend + 1;
    end
    t = t+dt;
    if floor(50*t)-floor(50*(t-dt))==1
        plot(iceField.x,iceField.y,ship.x,ship.y)
        xlim([-100 100]/scale)
        ylim([-50 50]/scale)
        title(num2str(t))
        vk(4)
        ri(3:4)
%         F_rot_acc = F_rot_acc+F_rotate(4)
        hold on
        if isfield(iceField,'icepiece')
            ind_rot_active =  iceField.icepiece.IndRot == 1;
            scatter(iceField.icepiece.loc(ind_rot_active,1),iceField.icepiece.loc(ind_rot_active,2),'MarkerEdgeColor','r','MarkerFaceColor','r')
            ind_rot_inactive =  iceField.icepiece.IndRot == 0;
            scatter(iceField.icepiece.loc(ind_rot_inactive,1),iceField.icepiece.loc(ind_rot_inactive,2),'MarkerEdgeColor','k','MarkerFaceColor','k')
        end
        pause(0.001)
        hold off
    end
    
    if abs(t - 76.65) < 0.0002
        disp('debug')
    end
       
    if floor(0.2*t)-floor(0.2*(t-dt))==1
        fprintf('t=%f\n',t)
    end
    %                     1  2     3     4     5     6      7    8      9            
    results(counter,:) = [t rk(1) rk(2) rk(3) rk(4) vk(1) vk(2) vk(3) vk(4) ...
        ak(1) ak(2) ak(3) ak(4) Fii(1) Fii(2) Fii(3) Fii(4) Rs(1) Ft(1) F_rotate(1) Fr(1) Rs(4) F_rotate(4)];
    %    10    11    12    13    14     15     16     17    18    19    20          21    22    23 
    counter = counter + 1;
%     disp(t)
end
   
    function [ iceForce, iceField, nLoad ] = iceForces( r, v, ship, iceField, vchan, t, lin,DB,areacorr)
        %This script determines ship ice contact and calls resolveContact to
        %calculate the component of the excitation force that is caused by ice
        %crushing and breaking
        
        %coordinate transformation: ice edge points to the body-fixed
        %coordinate system
        
        iceField.x = iceField.X*cos(r(4))+iceField.Y*sin(r(4))-(r(1)*cos(r(4))+r(2)*sin(r(4)));
        iceField.y = -iceField.X*sin(r(4))+iceField.Y*cos(r(4))+(r(1)*sin(r(4))-r(2)*cos(r(4)));
        xmid = [ship.mids,ship.mide];
        
%         plot(iceField.x,iceField.y,ship.x,ship.y)
%         axis([-100,100,-20,20])
        
        in = inhull(ship.x,ship.y,iceField.x,iceField.y,max(ship.x),min(ship.x),max(ship.y),min(ship.y),xmid);
        reallyin = find(in);
        %determine wether waterline and ice contact
        if isempty(find(reallyin, 1))
            iceForce = [0 0 0 0]';
            nLoad = [];
        else
            [iceForce, iceField, nLoad] = resolveContact(r, v, ship, iceField, reallyin, vchan, t, lin, DB, areacorr);
        end
        
    end

% This part has no use any more
%     function Fsub = submersion(deltarho,meanh,B,T,mu,L,tanphi,tanalpha,cosphi,sinphi,cospsi)
%         Fbow_buo = deltarho*9.81*meanh*B*T*(B+T)/(B+2*T);
%         Fbow_fric = deltarho*9.81*meanh*B*mu*(0.7*L - T/tanphi - B/(4*tanalpha));
%         Fmid = deltarho*9.81*meanh*B*mu*T*cosphi*cospsi*sqrt(1/sinphi^2+1/tanalpha^2);
%     end

    function [Frotate, iceField] = rotateForces(iceField, ship, r,v,rhoi,rhow,dt,mu,t, Hbw, xc, k)
        A_wave = (sqrt(1+2*k*Hbw)-1)/k;
        icepiece = iceField.icepiece;
        if icepiece.num ~= 0
            if floor(t/0.01)-floor((t+dt)/0.01) == -1
                dt = 0.01;
                del = nan(1,icepiece.num);
                alpha = nan(1,icepiece.num);
                phi = nan(1,icepiece.num);
                j = 0;
                for i = 1:icepiece.num
                    deleted = 0;
                    midx = icepiece.loc(i,1);
                    xi = ship.stem-midx;
                    itai = max(A_wave*cos(k*(xi-xc)) + 0.5*k*A_wave^2*cos(2*k*(xi-xc)),0);
                    if isinf(itai)
                        disp('itai is inf')
                    end
    %                 tic
                    ind = find(ship.x(length(ship.x)/2+1:end)<midx,1);
                    midy = -ship.y(length(ship.x)/2+ind)*sign(icepiece.loc(i,2));
                    icepiece.loc(i,2) = midy;
                    if icepiece.pos(i) == 1
                        alpha(i) = pi/2;
                    else
                        alpha(i) = atan(ship.tangenty(length(ship.x)/2+ind)/ship.tangentx(length(ship.y)/2+ind));
                    end
                    phi(i) = ship.phi(length(ship.x)/2+ind);
                    % midship
                    if midx>ship.mids && midx<ship.mide
                        phi(i) = phi(i) - abs(r(3));
                    end
    %                 phi(i) = atan(interp1(ship.x(length(ship.x)/2+1:end),ship.phi(length(ship.x)/2+1:end),icepiece.loc(i,1)));
                    dphi_over_dt = (1-icepiece.pos(i))*v(1)*ship.dphi_over_dx(length(ship.x)/2+ind);
    %                 dphi_over_dt = v(1)*interp1(ship.x(length(ship.x)/2+1:end),ship.dphi_over_dx(length(ship.x)/2+1:end),icepiece.loc(i,1));
                    d2phi_over_dt2 = (1-icepiece.pos(i))*v(1)*ship.d2phi_over_dx2(length(ship.x)/2+ind);
    %                 d2phi_over_dt2 = v(1)*interp1(ship.x(length(ship.x)/2+1:end),ship.d2phi_over_dx2(length(ship.x)/2+1:end),icepiece.loc(i,1));
                    vrel = [v(1)-v(4).*midy v(2)+v(4).*midx-v(3).*ship.zarm];
                    vreln = dot([vrel(1),vrel(2)*sign(icepiece.loc(i,2))], [sin(alpha(i)),cos(alpha(i))],2);
                    vrelt = dot([vrel(1),vrel(2)*sign(icepiece.loc(i,2))], [cos(alpha(i)),-sin(alpha(i))],2);
                    icepiece.RelDis(i) = min(0,icepiece.RelDis(i) + vreln*dt);
                    
                    
                    dvn_over_dt = 0;
                    L_i = icepiece.L(i);
                    S = icepiece.S(i);
                    Sc = icepiece.Sc(i);
                    SrhoBsq = icepiece.SrhoBsq(i);
                    hcons = icepiece.h(i)*(1-rhoi/rhow)/tan(phi(i));
                    mu_n(i) = mu*vreln*tan(phi(i))/sqrt(vrelt^2+(vreln*tan(phi(i)))^2);
                    mu_t(i) = mu*vrelt/sqrt(vrelt^2+(vreln*tan(phi(i)))^2);
                    mu_i = 0.05;
    %                 toc
    %                 tic
                    if icepiece.stage(i) == 0 && atan(icepiece.h(i)/L_i)+phi(i)>90/180*pi
                        icepiece.IndRot(i) = 0;
                    end
                    if icepiece.stage(i) == 0 && icepiece.IndRot(i) == 0 && atan(icepiece.h(i)/L_i)+phi(i)<90/180*pi && icepiece.RelDis(i) == 0
                        icepiece.IndRot(i) = 1;
                    end
                    if icepiece.IndRot(i) == 0
%                         j = j+1;
%                         del(j) = i; 
                        if icepiece.RelDis(i)==0
%                             Fn(i) = 0.5*icepiece.clength(i)*icepiece.h(i)*iceField.meansc;
                            Fn(i) = crushing(iceField.meansc, icepiece.clength(i)*icepiece.h(i),1);
                        else
                            Fn(i) = 0;
                        end
                    else
                        if icepiece.stage(i) == 0
        %                     tic
    %                         if atan(icepiece.h(i)/L_i)+phi(i)>85/180*pi
    %                             j = j+1;
    %                             del(j) = i;
    %                         else
                            [theta_initial,s,vx,vy,stage] = rotation_initialize_v2(vreln,dvn_over_dt,...
                                phi(i),dphi_over_dt,d2phi_over_dt2,rhoi,rhow,icepiece.h(i),L_i,dt,mu_n(i),mu_i,S,Sc,SrhoBsq,hcons,itai);
                            if isnan(s(1)) || isnan(theta_initial(1))
                                disp('theta is nan')
                            end
                            icepiece.theta(i,:) = theta_initial';
                            icepiece.s(i) = s;
                            icepiece.stage(i) = stage;
                            E_kinetic = 0.5*icepiece.m(3)*theta_initial(2)^2+0.5*icepiece.m(1)*(vx^2+vy^2);
                            E_potential = 0.25*rhow*9.81*Sc*Sc/S*(1-cos(2*theta_initial(1)));
                            E_initial = (E_kinetic+E_potential);%*(sin(phi(i))+mu_n*cos(phi(i)))^2;
                            icepiece.F_initial(i) = E_initial / (s*sin(phi(i))); % mean force during acceleration 
                            t_initial = s / vreln;
                            icepiece.t_initial(i) = t + t_initial;
                            if icepiece.F_initial(i)<0
                                pause(0.1)
                            end
    %                         end
        %                     toc
                        end
    %                     if atan(icepiece.h(i)/L_i)+phi(i)>85/180*pi
    %                         Fn(i) = 0;
    %                     else
                        [F(1:3),theta_next,x_next,y_next,s_next,stage] = rotation_v2(icepiece.m(i,:),icepiece.theta(i,:),...
                            vreln,dvn_over_dt,phi(i),dphi_over_dt,d2phi_over_dt2,rhoi,rhow,icepiece.h(i),L_i,dt,mu_n(i),mu_i,...
                            S,Sc,SrhoBsq,icepiece.y(i,:),icepiece.x(i,:),icepiece.s(i),hcons,icepiece.stage(i),itai);
                        if isnan(x_next(1)) || isnan(theta_next(1))
                            disp('theta is nan')
                        end
%                         if icepiece.stage(i) == 2
%                             pause(0.001)
%                         end
                        if F(1)>1e8
                            pause(0.1)
                        end
                        if F(1)<0
                            j = j+1;
                            del(j) = i;
                            deleted = 1;
                            F(1) = 0;
                            if icepiece.loc(i,1) < ship.mide && icepiece.loc(i,1) > ship.mids
                                disp('midship')
                            end
                        end
                        
                        icepiece.theta(i,:) = theta_next';
                        icepiece.x(i,:) = x_next';    
                        icepiece.y(i,:) = y_next'; 
                        icepiece.s(i) = s_next;
                        icepiece.stage(i) = stage;
                        icepiece.Fn(i) = F(1);
                        Fn(i) = F(1);
                        if t<icepiece.t_initial(i)
                            Fn(i) = Fn(i) + icepiece.F_initial(i);
                        end

                        if theta_next(1) > phi(i)
                            if deleted == 0
                                j = j+1;
                                del(j) = i;
                                deleted = 1;
                            end
                            vx_end = theta_next(2)*Sc/S*sin(phi(i))-icepiece.x(i,2);
                            vy_end = icepiece.y(i,2)-theta_next(2)*Sc/S*cos(phi(i));
                            vn_end = vx_end*sin(phi(i))-vy_end*cos(phi(i));
                            E_slamming = max(icepiece.m(i,1)*(vreln*sin(phi(i))-vn_end)*vreln*sin(phi(i)),0);
                            Fn(i) = Fn(i) + E_slamming/(vreln*dt*sin(phi(i)));
                        end
                    end
      
                    icepiece.loc(i,1) = icepiece.pos(i)*icepiece.loc(i,1) + ...
                            (1-icepiece.pos(i))*(icepiece.loc(i,1)-dt*v(1));
                    if icepiece.loc(i,1)<ship.min
                        if deleted == 0
                            j = j+1;
                            del(j) = i;
                            deleted = 1;
                        end
                    end
                end
                Fx = -Fn.*(sin(phi)+mu_n.*cos(phi)).*sin(alpha) - mu_t.*Fn.*cos(alpha);
                Fy = -sign(icepiece.loc(:,2)').*Fn.*(sin(phi)+mu_n.*cos(phi)).*cos(alpha) + sign(icepiece.loc(:,2)').*mu_t.*Fn.*sin(alpha); %CHECK SIGN
%                 Fz = Fn.*(cos(phi)-mu_n.*sin(phi));
%                 Mx = Fz.*icepiece.loc(:,2)'-Fy.*ship.zarm;
                Mx = -Fy.*ship.zarm;
%                 My = -Fz.*icepiece.loc(:,1)'+Fx.*ship.zarm;
                Mz = -Fx .* icepiece.loc(:,2)' + Fy .* icepiece.loc(:,1)';
                Frotate = [sum(Fx(isfinite(Fx))); sum(Fy(isfinite(Fy))); sum(Mx(isfinite(Mx))); sum(Mz(isfinite(Mz)))];
%                     [sum(Fx(isfinite(Fx))); sum(Fy(isfinite(Fy))); sum(Fz(isfinite(Fz))); sum(Mx(isfinite(Mx))); sum(My(isfinite(My))); sum(Mz(isfinite(Mz)))];
%                 Frotate = [sum(Fx(isfinite(Fx))); 0; 0;0; 0; 0];

                icepiece.Frotate = Frotate;
                %             toc
                del(isnan(del)) = [];
    %             if ~isempty(del)
    %                 icepiece.num
    %             end
                icepiece.num = icepiece.num-j;
%                 icepiece.num
                icepiece.RelDis(del) = [];
                icepiece.clength(del) = [];
                icepiece.IndRot(del) = [];
                icepiece.loc(del,:) = [];
                icepiece.L(del) = [];
                icepiece.Lprime(del) = [];
                icepiece.thetaw(del) = [];
                icepiece.S(del) = []; % 5
                icepiece.Sc(del) = [];
                icepiece.SrhoBsq(del) = [];
                icepiece.stage(del) = [];
                icepiece.vn(del) = [];
                icepiece.s(del) = []; % 10
                icepiece.m(del,:) = [];
                icepiece.theta(del,:) = [];
                icepiece.x(del,:) = [];
                icepiece.y(del,:) = [];
                icepiece.h(del) = [];
                icepiece.F_initial(del) = [];
                icepiece.Fn(del) = [];
                icepiece.t_initial(del) = [];
                icepiece.pos(del) = [];% 16
                iceField.icepiece = icepiece;
            else
                Frotate = icepiece.Frotate;
            end
        else
            Frotate = [0,0,0,0]';
        end
        if length(iceField.icepiece.loc(:,1)) ~= iceField.icepiece.num
            error('unequal number')
        end
    end

end
