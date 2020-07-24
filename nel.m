function Nel = nel(x,der)

for iii = 1:length(x)
    product = 0;
    for k = 1:20
        c = 1;
        cc = factorial(4*k+1)/factorial(4*k+1-der);
        for n = 1:k
            c = c*(4*n+1)*(4*n)^2*(4*n-1);
        end
        d = 0;
        for r = 1:k
            d = d+(1/(4*r-1)+0.5/r+1/(4*r+1));
        end
        cd = (-1)^k/c*d;
        product = product + cc*cd*x(iii)^(4*k+1-der);
    end

    if x(iii) == 0 && der == 0
        y2 = 0;
    elseif der == 0
        y2 = nev(x(iii),1,0)*log(x(iii));
    elseif der == 1
        y2 = nev(x(iii),1,1)*log(x(iii))+nev(x(iii),1,0)/x(iii);
    elseif der == 2
        y2 = nev(x(iii),1,2)*log(x(iii))+2*nev(x(iii),1,1)/x(iii)-nev(x(iii),1,0)/x(iii)^2;
    elseif der == 3
        y2 = nev(x(iii),1,3)*log(x(iii))+nev(x(iii),1,2)*3/x(iii)-nev(x(iii),1,1)*3/x(iii)^2+nev(x(iii),1,0)*2/x(iii)^3;
    end
    Nel(iii) = y2-product;
end