function [a,b] = solveEq(r, p, LoadCenter)

% lc = (E*h^3/12/1025/9.81)^0.25;

% No use any more
% load du3data.mat
% dd_Du3 = interp1(ddDu3(:,1),ddDu3(:,2),r);
% ddd_Du3 = interp1(dddDu3(:,1),dddDu3(:,2),r);

Coeff = [du2(r,2) du3(r,2); du2(r,2)+r*du2(r,3) du3(r,2)+r*du3(r,3)];
%Coeff = [du2(r,2) du3(r,2); r*du2(r,3) r*du3(r,3)];
Force = [LoadCenter*p*r*r/2/1025/9.81;p*r*r/2/1025/9.81];

result = inv(Coeff)*Force;

a = result(1);
b = result(2);