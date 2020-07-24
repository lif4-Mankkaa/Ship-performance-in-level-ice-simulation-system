function Nevm = nev(x,m,der)

for iii = 1:length(x)
    if m<der
        Nevm(iii) = 0;
    else
        Nevm(iii) = factorial(m)/factorial(m-der)*x(iii)^(m-der);
    end

    for k = 1:20
        Pi = 1;
        for n = 1:k
            Pi = Pi*(m+4*n)*(m+4*n-1)^2*(m+4*n-2);
        end
        coef = (-1)^k/Pi; 
        if m+4*k-der<0
            coef = 0;
        else
            coef = factorial(m+4*k)/factorial(m+4*k-der)*coef;
        end
        Nevm(iii) = Nevm(iii) + coef*x(iii)^(m+4*k-der);
    end
end