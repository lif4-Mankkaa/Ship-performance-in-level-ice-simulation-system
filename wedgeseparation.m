function [xy_plus, xy_minus] = wedgeseparation(x,y,filter)
% This is used as the 1st step to find 'b' in midship contact
% x is the x coordinate of ice field in ship coordinate system
% y is the y coordinate of ice field in ship coordinate system
% filter is the starting and end point of midship section

ind_filter = (x>=filter(1) & x<=filter(2));
x = x(ind_filter);
y = y(ind_filter);

indplus = (y>0);
indminus = (y<0);

xplus = x(indplus);
yplus = y(indplus);
xminus = x(indminus);
yminus = y(indminus);

j = 1;
k = 1;
ind_plus = [];
ind_minus = [];
for i = 2:length(xplus)-1
    if (yplus(i)>0) && (yplus(i)>yplus(i-1)) && (yplus(i)>yplus(i+1))
        ind_plus(j) = i;
        j = j+1;
    end
end
for i = 2:length(xminus)-1
    if (yminus(i)<0) && (yminus(i)<yminus(i-1)) && (yminus(i)<yminus(i+1))
        ind_minus(k) = i;
        k = k+1;
    end
end

x_plus = [filter(1),xplus(ind_plus),filter(2)];
y_plus = [mean(yplus(ind_plus)),yplus(ind_plus),mean(yplus(ind_plus))];
x_minus = [filter(1),xminus(ind_minus),filter(2)];
y_minus = [mean(yminus(ind_minus)),yminus(ind_minus),mean(yminus(ind_minus))];

xy_plus = [x_plus(:),y_plus(:)];
xy_minus = [x_minus(:),y_minus(:)];