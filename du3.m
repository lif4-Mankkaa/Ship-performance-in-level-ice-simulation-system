function Du3 = du3(x,der)

    Du3 = pi^1.5/4/gamma(0.75)^2*nev(x,0,der)-(log(exp(1)*sqrt(2))-0.577)*nev(x,1,der)-...
        gamma(0.75)^2/4/pi^0.5*nev(x,2,der)+nel(x,der);


end