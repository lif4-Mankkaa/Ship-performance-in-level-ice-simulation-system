function Du2 = du2(x,der)

% y(1) = 1;
% y(2) = -0.8472;
% y(3) = 0.2285;
% y(4) = 0;
% 
% for i = 5:100
%     y(i) = -y(i-4)/((i-1)*(i-2)*(i-2)*(i-3));
% end
% 
% Du2 = 0;
% for i = 1:100
%     if der>0
%         c = factorial(i-1+der)/factorial(i-1);
%     else
%         c = 1;
%     end
%     if i + der <= 100
%         Du2 = Du2 + y(i+der)*x^(i-1)*c;
%     end
% end
% 
% k = Du2 * pi^1.5/2/gamma(0.75)^2
Du2 = nev(x,0,der)*pi^1.5/2/gamma(0.75)^2 - pi/2*nev(x,1,der) + gamma(0.75)^2/2/pi^0.5*nev(x,2,der);