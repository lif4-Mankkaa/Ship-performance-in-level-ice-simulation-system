function keix = kei(x)
%https://en.wikipedia.org/wiki/Kelvin_functions
term1 = -log(x/2)*bei(x);
term2 = -pi/4*ber(x);

term3 = nan(1,10);
for k = 0:1:9
    term3(k+1) = (-1)^k*psi(2*k+2)/factorial(2*k+1)*(x^2/4)^(2*k+1);
end
conv = cumsum(term3);
term3 = sum(term3);

keix = term1 + term2 + term3;

    function beix = bei(x)
        
    end
end