function sigma0 = stressbow(h,r,der)
E = 5e9;
rho = 1025;
g = 9.81;
nu = 0.33;

sigmac = 1.28e6;
phi = 35;
q = sigmac*cosd(phi);
q = 10;

k = rho*g;
D = E*h^3/12/(1-nu^2);
p = q*pi*r^2/2;
x0 = 4*r/3/pi;
beta = (k/D)^0.25;
l = 1/beta;

% kei = imag(besselk(0,(r/l)*(1+1i)/sqrt(2)));
% kei_delta1 = imag(besselk(0,(r/l+r/l/1000)*(1+1i)/sqrt(2)));
% kei_delta2 = imag(besselk(0,(r/l-r/l/1000)*(1+1i)/sqrt(2)));
kei1 = interp1(der(:,1),der(:,2),r/l);

term1 = 3*q*r*D*(1+nu)*kei1/(k*h^2*l^3);

interval = 0.1;
xi = 0:interval:2E4;
Q = sqrt( 0.5*sqrt(xi.^4+beta.^4)+xi.^2/2 );
% K = sqrt( max(0.5*sqrt(xi.^4+beta.^4)-xi.^2/2,0) );
K = 1/sqrt(2)*beta^2./sqrt((sqrt(xi.^4+beta^4)+xi.^2));
Qstar = Q.*(2*K.^2+(1-nu)*xi.^2)./(K.*(2*Q.^2-(1-nu)*xi.^2));
KoverQ = beta^2./((xi.^4+beta^4).^0.5+xi.^2);
Qstar = (beta^2*(KoverQ)+(1-nu)*xi.^2) ./ (beta^2-(1-nu)*xi.^2.*(KoverQ));
K2plusQ2 = (xi.^4+beta^4).^0.5;
A1 = p*beta^2/pi/k*...
            exp(-Q*x0)./sqrt(xi.^4+beta^4).*...
            (2*Q.^2-(1-nu)*xi.^4)./(4*Q.^2.*(K.^2+(1-nu)*xi.^2)-(1-nu)^2*xi.^4).*...
            ( K.*(K2plusQ2+nu*xi.^2).*cos(K*x0)-Q.*(K2plusQ2-nu*xi.^2).*sin(K*x0) );

%Q*K = 0.5beta^2
term2_int = A1.*(nu*(xi.^2+beta^2.*Qstar)-xi.^2);
% plot(cumsum(term2_int)); 
term2 = interval*sum(term2_int)*6*D/h^2;
sigma0 = term1 + term2;