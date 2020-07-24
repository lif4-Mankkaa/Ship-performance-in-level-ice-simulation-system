function Fr = rudderforce(Rasp, alpha, chord, speed, Ar, rudpos, wake, Kr, T, D)
%This is to calculate rudder force
%Rasp: aspect ratio
%alpha: attack angle
%chord: chord length of rudder
%speed: water speed coming to rudder
%Ar: rudder area
option = 1;
nu = 1.35e-6;
rho = 1025;
if option == 1
    v0 = speed*(1-wake);
    Kt_over_J2 = T/(v0^2*1025*D^2);
    vr = v0 + Kr*v0*(sqrt(1+8*Kt_over_J2/pi)-1);
    vr = 0.80*vr;
    % vr = 1.2*0.85*speed;

    Rasp = 1.3*Rasp;

    Cq = 1; %if no better determination
    Cl = 2*pi * Rasp*(Rasp+1)/(Rasp+2)^2 * sind(alpha);
    % Cl = 2*pi * Rasp*(Rasp+0.7)/(Rasp+1.7)^2 * sind(alpha) + Cq*sind(alpha)...
    %     * abs(sind(alpha)) * cosd(alpha);

    Rn = vr*chord/nu;

    Cd0 = 2.5*0.075/(log(Rn)-2)^2;

    Cd = 1.1*Cl^2/(pi*Rasp)+ Cd0;% + Cq*abs(sind(alpha))^3 + Cd0;

    Fr = [0;0;0;0];

    Fr(1) = -Cd*0.5*rho*vr^2*Ar;
    Fr(2) = Cl*0.5*rho*vr^2*Ar;
    Fr(3) = -Fr(2)*rudpos(2);
    % Fr(5) = -Fr(1)*rudpos(2);
    Fr(4) = rudpos(1)*Fr(2);
elseif option == 2
    Cl = 2*pi*Rasp*(Rasp+0.7)/(Rasp+1.7)^2*sind(alpha) + sind(alpha)*abs(sind(alpha))*cosd(alpha);
    v0 = speed*(1-wake);
    Kt_over_J2 = T/(v0^2*1025*D^2);
    deltaL = T*sind(alpha)*(1+1/sqrt(1+8*Kt_over_J2/pi));
    Rn = v0*chord/nu;
    Cd0 = 2.5*0.075/(log(Rn)-2)^2;
    Cd = 1.1*Cl^2/(pi*Rasp)+ Cd0;% + Cq*abs(sind(alpha))^3 + Cd0;
    deltaD = T*(1-cosd(alpha))*(1+1/sqrt(1+8*Kt_over_J2/pi));
    Fr = [0;0;0;0];
    Fr(1) = -Cd*0.5*rho*v0^2*Ar-deltaD;
    Fr(2) = Cl*0.5*rho*v0^2*Ar+deltaL;
    Fr(3) = -Fr(2)*rudpos(2);
    % Fr(5) = -Fr(1)*rudpos(2);
    Fr(4) = rudpos(1)*Fr(2);
end
Fr = 1.2*Fr;