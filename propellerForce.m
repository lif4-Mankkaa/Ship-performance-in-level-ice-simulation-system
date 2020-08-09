function [ propForce ] = propellerForce( v, ship, pp_D, power_t, prop_type )
%calculation of propeller forces
%pp_D: propeller pitch to diameter ratio, if fixed pitch, pp_D=0
%power_t: ship power to propeller at time t
%prop_type: blade area ratio, 70, 85 or 100
K_ref = -0.5848*ship.refpitch+3.349;
Ke_ref = (0.250/pi^2)^(1/3)*K_ref;
refbollard = 2*Ke_ref*(1/2*ship.refpower*ship.dprop)^(2/3)*1e3;

if prop_type == 70
    C_ita = [-0.3742 1.054 -0.003];
    C_ke = [-0.6512 3.3928];
elseif prop_type == 85
    C_ita = [-0.3304 0.9582 0.0176];
    C_ke = [-0.5848 3.349];
elseif prop_type == 100
    C_ita = [-0.3839 1.0854 -0.0827];
    C_ke = [-0.3033 3.0321];
else
    error('prop_type should be 70 85 or 100')
end

itaref = C_ita(1)*(ship.refpitch)^2+C_ita(2)*ship.refpitch +C_ita(3);

propForce = [0 0 0 0]';
powerratio = power_t/ship.fullpower;
ita = C_ita(1)*(pp_D)^2+C_ita(2)*pp_D +C_ita(3);
itaratio = (pp_D == 0)*1 + (pp_D ~= 0)*ita/itaref;

vow = (powerratio*itaratio)^0.3064*ship.vow;
if pp_D ~= 0
    bollard = powerratio^0.6666*(C_ke(1)*pp_D+C_ke(2))/(C_ke(1)*ship.refpitch+...
        C_ke(2))*refbollard;
else
    bollard = powerratio^0.6666*refbollard;
end

if 0 <= v(1) && v(1) < vow
    propForce(1) = bollard * (1-v(1)/(3*vow)-2/3*(v(1)/vow)^2);
elseif v(1) < 0
    propForce(1) = bollard;
end

end

