for jj = 1:5000
    ii(jj) = imag(besselk(0,(jj*0.0001)*(1+1i)/sqrt(2)));
end

der = (ii(3:end) - ii(1:end-2))/0.0002;
der = [[0.0002:0.0001:0.4999]',der(:)];



ss = [];

for mm = 1:50
    ss(mm) = stressbow(0.3,mm*0.01,der);
end

figure
plot(0.01:0.01:0.5,ss)