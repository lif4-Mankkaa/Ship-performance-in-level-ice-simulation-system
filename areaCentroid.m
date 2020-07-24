function [ midx, midy ] = areaCentroid( x0, y0, theta, icex, icey, wlx, wly )
%centroid of an area bounded by polylines and bisected by a span line.
x0 = x0(:);
y0 = y0(:);
wlx = wlx(:);
wly = wly(:);

coeff = 1/(cos(theta)^2 + sin(theta)^2);
icet = coeff * (icex * cos(theta) + icey * sin(theta)) - (x0 * cos(theta) + y0 * sin(theta));
icen = coeff * (-icex * sin(theta) + icey * cos(theta)) + (x0 * sin(theta) - y0 * cos(theta));
wlt = coeff * (wlx *  cos(theta) + wly * sin(theta)) - (x0 * cos(theta) + y0 * sin(theta));
wln = coeff * (-wlx * sin(theta) + wly * cos(theta)) + (x0 * sin(theta) - y0 * cos(theta));

%centroid of the area bounded by the span line and waterline
diffx = diff(icet);
for i = 1:length(icet)-1
    sumy1(i) = 2*icen(i+1)+icen(i);
    sumy2(i) = icen(i)+icen(i+1);
end
Si = 0.5*sumy2.*abs(diffx');
if sumy2 == 0
    meant_i = mean(wlt);
    meann_i = mean(wln);
else
    sumy2(sumy2==0) = Inf;
    Cix = 1/3*diffx'.*sumy1./sumy2+icet(1:end-1)';
    if sum(abs(Si)) == 0
        meant_i = 0;
        meann_i = 0;
    else
        meant_i = sum(abs(Si).*Cix)/sum(abs(Si));
        Ciy = (0.5*min(abs(icen(1:end-1)'),abs(icen(2:end)')).^2+...
            (min(abs(icen(1:end-1)'),abs(icen(2:end)'))+1/3*abs(diff(icen')))*0.5.*abs(diff(icen')))./...
            (min(abs(icen(1:end-1)'),abs(icen(2:end)'))+0.5.*abs(diff(icen')));
        meann_i = sum(Si.*Ciy)/sum(abs(Si));
    end
end

diffx = diff(wlt);
sumy1 = [];
sumy2 = [];
for i = 1:length(wlt)-1
    sumy1(i) = 2*wln(i+1)+wln(i);
    sumy2(i) = wln(i)+wln(i+1);
end
Si = 0.5*sumy2.*abs(diffx');
if sumy2 == 0
    meant_w = mean(wlt);
    meann_w = mean(wln);
else
    sumy2(sumy2==0) = Inf;
    Cix = 1/3*diffx'.*sumy1./sumy2+wlt(1:end-1)';
    if sum(abs(Si)) == 0
        meant_w = 0;
        meann_w = 0;
    else
        meant_w = sum(abs(Si).*Cix)/sum(abs(Si));
        Ciy = (0.5*min(abs(wln(1:end-1)'),abs(wln(2:end)')).^2+...
            (min(abs(wln(1:end-1)'),abs(wln(2:end)'))+1/3*abs(diff(wln')))*0.5.*abs(diff(wln')))./...
            (min(abs(wln(1:end-1)'),abs(wln(2:end)'))+0.5.*abs(diff(wln')));
        meann_w = sum(Si.*Ciy)/sum(abs(Si));
    end
end
% meant_w = mean(vertcat(wlt, wlt(2:end-1)));
% meann_iw = zeros(length(wln)-3, 1);
% 
% for i = 2:length(wln)-2
%     meann_iw(i-1) = mean([wln(i:i+1); 0; 0]);
% end
% 
% meann_sw = mean([0; wln(2); 0]);
% meann_ew = mean([0; wln(end-1); 0]);
% meann_w = mean([meann_sw; meann_iw; meann_ew]);

%centroid of the area bounded by span line and ice edge. 
% meant_i = mean(vertcat(icet, icet(2:end-1)));
% meanni_i = zeros(length(icen)-3, 1);

% for i = 2:length(icen)-2
%     meanni_i(i-1) = mean([icen(i:i+1); 0; 0]);
% end
% 
% meann_si = mean([0; icen(2); 0]);
% meann_ei = mean([0; icen(end-1); 0]);
% meann_i = mean([meann_si; meanni_i; meann_ei]);

%centroid of whole area
S_ice = trapz(icet,icen);
S_wl  = trapz(wlt,wln);
meant = (meant_w*S_wl+meant_i*S_ice)/(S_wl+S_ice);
if S_wl+S_ice == 0
    meann = mean([meann_w,meann_i]);
    meant = mean([meant_w,meant_i]);
else
    if sign(icen(2:end-1)) == sign(wln(2))
        meann = (S_ice>S_wl)*(S_ice*meann_i - S_wl*meann_w)/(S_wl+S_ice)+...
            (S_ice<=S_wl)*(S_wl*meann_w - S_ice*meann_i)/(S_wl+S_ice);
    else
        meann = (S_ice*meann_i + S_wl*meann_w)/(S_wl+S_ice); 
    end
end

%coordinate transformation
midx = x0 + meant*cos(theta)-meann*sin(theta);
midy = y0 + meant*sin(theta)+meann*cos(theta);

end

