function [Corr_v, Corr_A, Corr_ind] = areacorrec(ind0,h,p,psi,E)
% This function aims to find the correction factore gamma to modify the
% contact area
% h: ice thickness, m
% p: nominal pressure, Pa
% psi: flare angle
% ind0: indentation depth
% E: elastic modulus, Pa
% tic
lc = ((E*h^3)/12/0.91/1025/9.81)^0.25;
%ind_init = ind0;
def0 = 0;
e = 1;
def_inc = ind0*tand(psi); %deflection increment
dtheta = 1;
iter = 0;
while e>0.01
    % Calculate area before correction
    ind = ind0 - def0/tand(psi);
    if h/tand(psi) > ind
        LC = 1/3;
        A0 = 0.5*ind*ind*dtheta;
    else
        a_trap = ind-h/tand(psi);
        b_trap = ind;
        LC = (h/tand(psi)/3*(2*a_trap+b_trap)/(a_trap+b_trap))/ind;
        A0 = 0.5*(a_trap+b_trap)*dtheta*h/tand(psi);
    end

    r = ind/lc;
    [a,b] = solveEq(r, p, LC);
    def = a*du2(0,0)+b*du3(0,0);
    
%     if def0 == 0
%         %record def
%         def1 = def;
%     end
    
    e = abs(def-def0)/def0;
    if e >0.01
        def_inc = 0.5*def_inc;
        if def > def0
            def0 = def0 + def_inc;
        else
            def0 = def0 - def_inc;
        end
    end
    iter = iter+1;
    if iter>1000
        error('Too many iteration in area correction')
    end
end

Corr_v = def0/(ind0*tand(psi));
A_corr = A0;
if h/tand(psi) > ind0
    A0 = 0.5*ind0*ind0*dtheta;
else
    a_trap = ind0-h/tand(psi);
    b_trap = ind0;
    A0 = 0.5*(a_trap+b_trap)*dtheta*h/tand(psi);
end
Corr_A = A_corr/A0;
Corr_ind = ind/ind0;

% toc