function Fcr = crushing(sigmac, area, lin)       

if lin == 0
    Fcr = sigmac.*(area/0.05^2).^(-0.3).*area;
elseif lin == 1
    %linear
    Fcr = sigmac.*area;
end