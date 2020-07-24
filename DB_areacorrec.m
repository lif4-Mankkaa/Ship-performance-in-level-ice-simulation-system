function [V,A,ind] = DB_areacorrec(mean_h,thrh_h,mean_sigc,thrh_sigc,mean_E,thrh_E)

lc = ((mean_E*mean_h^3)/12/0.91/1025/9.81)^0.25;
r = 0.001:0.01:0.201;
depth = r*lc;
if thrh_h ~= 0
    h = linspace(mean_h-thrh_h,mean_h+thrh_h,11);
else
    h = mean_h;
end
if thrh_sigc ~= 0
    sigc = linspace(mean_sigc-thrh_sigc,mean_sigc+thrh_sigc,11);
else
    sigc = mean_sigc;
end
if thrh_E ~= 0 
    E = linspace(mean_E-thrh_E,mean_E+thrh_E,11);
else
    E = mean_E;
end
meanPhi = 20:5:80;
for i = 1:length(r)
    for j = 1:length(h)
        for k = 1:length(sigc)
            for m = 1:length(meanPhi)
                for n = 1:length(E) 
                    [Corr_v, Corr_A, Corr_ind] = areacorrec(depth(i),h(j),sigc(k),meanPhi(m),E(n));
                    V.depth(i).h(j).sigc(k).phi(m).E(n) = Corr_v;
                    A.depth(i).h(j).sigc(k).phi(m).E(n) = Corr_A;
                    ind.depth(i).h(j).sigc(k).phi(m).E(n) = Corr_ind;
                end
            end
        end
    end
end
    