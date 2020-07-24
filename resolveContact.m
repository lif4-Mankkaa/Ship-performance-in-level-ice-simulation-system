function [ iceForce, iceField, nLoad ] = resolveContact( r, v, ship, iceField, reallyin, vchan, ~, lin, DB, areacorr)
%resolves iceContacts up to bending failure and alters the ice field edge if bending failure occurs
%/////////////////////////////////////////////////////////////////////////
%vector reallyin containing the indices of ice nodes inside (not merely on)
%the waterline polygon is manipulated to ensure the grouping algorithm
%works. the grouping algorithm compiles the vector contacts with
%length
%twice the number of discrete ice contacts. 'contacts' contains as odd
%indiced values the index of the first ice node inside the contact and as
%even indiced values the index of the last ice node inside the contact. if
%the contact is comprised of only one ice node, the same index is given as
%the first and the last in contact.
cy3 = cos(r(4));
sy3 = sin(r(4));

[contacts, noOfContacts] = grouping(reallyin);
% noOfContacts
%% /////////////////////////////////////////////////////////////////////////
%/////////////////////////////////////////////////////////////////////////
%determination of intersection points for each contact
ice_debug1 = iceField;
%outside contains indices of the ice nodes just outside (or on) the hull
%polygon
outside = contacts + repmat([-1 1], 1, noOfContacts);

%go through each contact and determine intersection points and hull
%nodes in the contact
contactswl = zeros(1,noOfContacts*2);
contactswl_1 = zeros(1,noOfContacts*2);
xint = zeros(noOfContacts*2,1);
yint = zeros(noOfContacts*2,1);
midx = zeros(noOfContacts,1);
midy = zeros(noOfContacts,1);
%         midX = zeros(noOfContacts,1);
%         midY = zeros(noOfContacts,1);

for i = 1:noOfContacts
    
    icexs = [iceField.x(outside(2*i-1)) iceField.x(contacts(2*i-1))];
    icexe = [iceField.x(outside(2*i)) iceField.x(contacts(2*i))];
    iceys = [iceField.y(outside(2*i-1)) iceField.y(contacts(2*i-1))];
    iceye = [iceField.y(outside(2*i)) iceField.y(contacts(2*i))];
    
    [xints, yints, iints] = polyxpoly(icexs, iceys, ship.x, ship.y);
    [xinte, yinte, iinte] = polyxpoly(icexe, iceye, ship.x, ship.y);
    
    if size(xints)>[1,0] & size(xinte) == [1,1]
        xints(xints == xinte) = [];
    elseif size(xinte)>[1,0] & size(xints) == [1,1]
        xinte(xinte == xints) = [];
    end
    
    %these ifs in case point is too close to the line for
    %interpolation to work
    
    if size(iinte) == [0,2]
        icexe = [iceField.x(outside(2*i)+1) iceField.x(contacts(2*i))];
        iceye = [iceField.y(outside(2*i)+1) iceField.y(contacts(2*i))];
        [xinte, yinte, iinte] = polyxpoly(icexe, iceye, ship.x, ship.y);
    end
    
    if size(iints) == [0,2]
        icexs = [iceField.x(outside(2*i-1)-1) iceField.x(contacts(2*i))];
        iceys = [iceField.y(outside(2*i-1)-1) iceField.y(contacts(2*i))];
        [xints, yints, iints] = polyxpoly(icexs, iceys, ship.x, ship.y);
    end
    
    if size(iints) == [0,2] | size(iinte) == [0,2]
        %ota kontakti pois
        noOfContacts = noOfContacts-1;
        xinte = NaN;
    end
    
    %in case only one intersection point
    if round(xinte,3) == round(xints,3) & round(yinte,3) == round(yints,3)
        xinte = NaN;
        noOfContacts = noOfContacts - 1;
    end
    
    if ~isnan(xinte)
        if iints(1,2)==iinte(1,2)
            contactswl(2*i-1:2*i) = [-10000 -10000];%no ship nodes inside
        else
            contactswl(2*i-1:2*i) = [iints(1,2) iinte(1,2)+1];
        end
        %in case no ship nodes inside contact
        contactswl_1(2*i-1:2*i) = [iints(1,2) iinte(1,2)+1];
        xint(2*i-1:2*i) = [xints(1) xinte(1)];
        yint(2*i-1:2*i) = [yints(1) yinte(1)];
        
        midx(i) = mean(xint(2*i-1:2*i));
        midy(i) = mean(yint(2*i-1:2*i));
    else
        contactswl(2*i-1:2*i) = [NaN NaN];
        contacts(2*i-1:2*i) = [NaN NaN];
        xint(2*i-1:2*i) = [NaN NaN];
        yint(2*i-1:2*i) = [NaN NaN];
        
        midx(i) = NaN;
        midy(i) = NaN;
    end
    %             midX(i) = y(1) + midx(i).*cy3 - midy(i).*sy3;
    %             midY(i) = y(2) + midx(i).*sy3 + midy(i).*cy3;
end

xint = round(xint,3);
yint = round(yint,3);

contactswl_1 = contactswl_1(~isnan(contactswl));
contactswl = contactswl(~isnan(contactswl));
contacts = contacts(~isnan(contacts));

xint = xint(~isnan(xint));
yint = yint(~isnan(yint));

midx = midx(~isnan(midx));
midy = midy(~isnan(midy));


%         if isempty(midx)
%             iceforce = 0;
%             iceField = iceField;
%             nload = [];
%         else
%% /////////////////////////////////////////////////////////////////////////
%/////////////////////////////////////////////////////////////////////////
%span line (tangent) equations of the form ax+by+c=0 for each contact are determined

%xdiff and ydiff contain the dx and dy of each contact.
xdiff = diff(xint);
xdiff = xdiff(1:2:end);
ydiff = diff(yint);
ydiff = ydiff(1:2:end);
cLength = sqrt(xdiff.^2+ydiff.^2);
%a, b and c are vectors of lenght n containing the coefficients for span
%lines of the form ax+by+c=0
a = ydiff./xdiff;
b = -1*ones(length(a),1);
c = -1*a.*xint(1:2:end-1)+yint(1:2:end-1);

a0 = ones(size(a));
if ~isempty(a(isinf(a)))
    b(isinf(a)) = 0;
    c(isinf(a)) = xint(2*find(isinf(a)));
    a(isinf(a)) = -1;
    a0(isinf(a)) = -inf;
end

%% /////////////////////////////////////////////////////////////////////////
%/////////////////////////////////////////////////////////////////////////
%interpolation of the properties of the ice field for the contacts
midX = r(1) + midx.*cy3 - midy.*sy3;
midY = r(2) + midx.*sy3 + midy.*cy3;
%ice thickness at midpoints is interpolated from the iceField
hi = interp1(iceField.N(100,:)-1, iceField.h(100,:), min(round(max(midX,0)),length(iceField.N(100,:))-1), 'nearest');
sigmaf = interp1(iceField.N(100,:)-1, iceField.sigmaf(100,:), min(round(max(midX,0)),length(iceField.N(100,:))-1), 'nearest');
sigmac = interp1(iceField.N(100,:)-1,iceField.sigmac(100,:), min(round(max(midX,0)),length(iceField.N(100,:))-1),'nearest');
l_c = (iceField.E*hi.^3/(12*(1-0.3.^2)*1025*9.81)).^0.25;

%% /////////////////////////////////////////////////////////////////////////
%/////////////////////////////////////////////////////////////////////////
%unit tangents and normals for each contact
tangent = zeros(noOfContacts,2);
normal = zeros(noOfContacts,2);

for i = 1:noOfContacts
    index = sort([contactswl_1(2*i-1) contactswl_1(2*i)]);
    dxt = sum(ship.tangentx(index(1):index(2)))./norm([sum(ship.tangentx(index(1):index(2))) sum(ship.tangenty(index(1):index(2)))]);
    dyt = sum(ship.tangenty(index(1):index(2)))./norm([sum(ship.tangentx(index(1):index(2))) sum(ship.tangenty(index(1):index(2)))]);
    dxn = sum(ship.normalx(index(1):index(2)))./norm([sum(ship.normalx(index(1):index(2))) sum(ship.normaly(index(1):index(2)))]);
    dyn = sum(ship.normaly(index(1):index(2)))./norm([sum(ship.normalx(index(1):index(2))) sum(ship.normaly(index(1):index(2)))]);
    tangent(i,:) = [dxt dyt];
    normal(i,:) = [dxn dyn];
end

meanPhi = zeros(noOfContacts,1);
for j = 1:noOfContacts
    %mean value of flare angle for the contact
    if contactswl_1(j*2-1) > contactswl_1(j*2)
        meanPhi(j) = mean(ship.phi(contactswl_1(2*j):contactswl_1(2*j-1)));
    else
        meanPhi(j) = mean(ship.phi(contactswl_1(2*j-1):contactswl_1(2*j)));
    end
    if midx(j)>ship.mids && midx(j)<ship.mide
        meanPhi(j) = meanPhi(j)-abs(r(3));
        %                r(4)
    end
end

%relative motion of the contact area midpoints. y(4:6) are the x-, y- and
%rotational component of vcog. vrelt, vreln and creln1 correspond to the
%like named forces in the paper.
vrel = [v(1)-v(4).*midy(:) v(2)+v(4).*midx(:)-v(3).*ship.zarm];
indforv = midy>0; % when channel is closing
vrel = vrel + [zeros(length(vrel(:,1)),1) indforv*vchan(1)-(1-indforv)*vchan(2)];
vrelt = dot(vrel(:,1:2), tangent,2);
vreln = dot(vrel(:,1:2), normal,2);
vreln1 = vreln .* cos(meanPhi);%+vrel(:,3).*sin(meanPhi);
vreln2 = vreln .* sin(meanPhi);%+vrel(:,3).*cos(meanPhi);

%% Determination of contact areas
%/////////////////////////////////////////////////////////////////////////
%/////////////////////////////////////////////////////////////////////////
%determination of contact areas

%initialization of area and meanPhi vectors
area = NaN(noOfContacts,1);
depth = NaN(noOfContacts,1);

cmPhi = cos(meanPhi);
smPhi = sin(meanPhi);

j = 1;

LC_DB = -1*ones(noOfContacts,1);
b_wedge = NaN(noOfContacts,1);
c_wedge = NaN(noOfContacts,1);
the = NaN(noOfContacts,1);

k_DB = [];
for i = 1:2:2*noOfContacts
    %waterline nodes inside contact
    %if a contact is across the bow
    
    if contactswl(i) == -10000 %no ship nodes inside
        xwl = [];
        ywl = [];
    elseif sign(ship.y(contactswl(i))) ~=sign(ship.y(contactswl(i+1)))
        if contactswl(i) > contactswl(i+1)
            xwl = ship.x(contactswl(i+1):1:contactswl(i));
            ywl = ship.y(contactswl(i+1):1:contactswl(i));
        elseif contactswl(i) < contactswl(i+1)
            xwl = vertcat(ship.x(contactswl(i+1):end), ship.x(1:contactswl(i)));
            ywl = vertcat(ship.y(contactswl(i+1):end), ship.y(1:contactswl(i)));
        end
    else
        if contactswl(i) > contactswl(i+1)
            xwl = ship.x(contactswl(i):-1:contactswl(i+1));
            ywl = ship.y(contactswl(i):-1:contactswl(i+1));
        elseif contactswl(i) < contactswl(i+1)
            xwl = ship.x(contactswl(i):contactswl(i+1));
            ywl = ship.y(contactswl(i):contactswl(i+1));
        else
            xwl = ship.x(contactswl(i));
            ywl = ship.y(contactswl(i));
        end
    end
    %ice nodes inside contact
    xice0 = iceField.x(contacts(i):contacts(i+1));
    yice0 = iceField.y(contacts(i):contacts(i+1));
    %distances from span line to ice points inside the contact area
    if ~isinf(a0(j))
        distance = (abs(a(j)*xice0+b(j)*yice0+c(j))/sqrt(a(j)^2+b(j)^2));
    else
        distance =(isinf(a(j))) .* (xice - midx(j));
    end
    
    correctedDistance = (distance.*abs(tan(meanPhi(j)))>hi(j)).*(hi(j)/sin(meanPhi(j))) + ...
        (distance.*abs(tan(meanPhi(j)))<=hi(j)).*distance./abs(cos(meanPhi(j)));
    %distance to move ice nodes for area calculation
    moveDistance = correctedDistance-distance;
    %determination of corrected ice nodes for area calculation
    
    if sign(ship.y(contactswl_1(i))) ==sign(ship.y(contactswl_1(i+1)))
        pos(j) = 0; % 0 for ship side
        if sign(a(j))*(a(j)*xice0+b(j)*yice0+c(j)) == 0
            xice = xice0 + moveDistance .* normal(j,1);
            yice = yice0 + moveDistance .* normal(j,2);
        else
            xice = logical(sign(a(j)).*(a(j)*xice0+b(j)*yice0+c(j)) > 0) .* (xice0 - moveDistance .* normal(j,1)) + ...
                logical(sign(a(j)).*(a(j)*xice0+b(j)*yice0+c(j)) < 0) .* (xice0 + moveDistance .* normal(j,1));
            yice = logical(sign(a(j)).*(a(j)*xice0+b(j)*yice0+c(j)) > 0) .* (yice0 - moveDistance .* normal(j,2)) + ...
                logical(sign(a(j)).*(a(j)*xice0+b(j)*yice0+c(j)) < 0) .* (yice0 + moveDistance .* normal(j,2));
        end
    else
        pos(j) = 1; % 1 for stem
        
        xice = xice0;
        yice = yice0;
    end
    
    if xice == 0 & yice == 0
        error('wo cao')%error
    end
    
    xice = xice(:);
    yice = yice(:);
    xwl = xwl(:);
    ywl = ywl(:);
    %collection of points to determine the contact area
    %             if sign(ship.y(contactswl(i))) ~=sign(ship.y(contactswl(i+1)))
    xbound = vertcat(xint(i), xice, xint(i+1), xwl, xint(i));
    ybound = vertcat(yint(i), yice, yint(i+1), ywl, yint(i));
    %             elseif midy(j) < 0
    %                 xbound = vertcat(xint(i), xice, xint(i+1), xwl, xint(i));
    %                 ybound = vertcat(yint(i), yice, yint(i+1), ywl, yint(i));
    %             elseif midy(j) >= 0
    %                 xbound = vertcat(xint(i), xice, xint(i+1), xwl, xint(i));
    %                 ybound = vertcat(yint(i), yice, yint(i+1), ywl, yint(i));
    %             end
    %refined midpoint of contact for possible alteration of ice field
    xin1 = vertcat(xint(i), xice, xint(i+1));
    yin1 = vertcat(yint(i), yice, yint(i+1));
    if contactswl(i) >= contactswl(i+1)
        xin2 = vertcat(xint(i+1), xwl(end:-1:1), xint(i));
        yin2 = vertcat(yint(i+1), ywl(end:-1:1), yint(i));
    elseif contactswl(i) < contactswl(i+1)
        xin2 = vertcat(xint(i+1), xwl, xint(i));
        yin2 = vertcat(yint(i+1), ywl, yint(i));
    end
    [ctrx(j), ctry(j)] = areaCentroid(xint(i), yint(i), atan(a(j)*a0(j)), xin1, yin1, xin2, yin2);
    if ctrx(j)<min(xint(i:i+1)) || ctrx(j)>max(xint(i:i+1)) || isnan(ctrx(j))
        ctrx(j) = mean(xint(i:i+1));
        ctry(j) = mean(yint(i:i+1));
    end
    %-------------------------------------------------------
    %crushing center and slope, crushing depth
    if meanPhi(j)/pi*180 <= 80
        center_dist = (abs(a(j)*ctrx(j)+b(j)*ctry(j)+c(j))/sqrt(a(j)^2+b(j)^2));%distance from center to crushing line
        c_DB = center_dist*cos(meanPhi(j));
        b_DB = max(distance);
        a_DB = min(hi(j)/tan(meanPhi(j)),b_DB);
        LC_DB(j) = (6*c_DB-3*a_DB)*c_DB/a_DB/(3*c_DB-2*a_DB);
        %             alpha1 = atand(abs((yint(i)-yice0(1))/(xint(i)-xice0(1))));
        %             alpha2 = atand(abs((yint(i+1)-yice0(end))/(xint(i+1)-xice0(end))));
        %             alpha = atand(abs(a(j)/b(j)));
        %
        %             alpha_unknown = -0.5*abs(alpha1-alpha2)+alpha;
        %             k_DB(j) = tand(alpha_unknown)

        %             alpha1 = atand((yint(i)-yice0(1))/(xint(i)-xice0(1)));
        %             alpha2 = atand((yint(i+1)-yice0(end))/(xint(i+1)-xice0(end)));
        %             alpha = atand(-a(j)/b(j));
    %     k1 = (yint(i)-yice0(1))/(xint(i)-xice0(1));
    %     k2 = (yint(i+1)-yice0(end))/(xint(i+1)-xice0(end));
        k1 = (iceField.y(contacts(i)-1)-yice0(1))/(iceField.x(contacts(i)-1)-xice0(1));
        k2 = (iceField.y(contacts(i+1)+1)-yice0(end))/(iceField.x(contacts(i+1)+1)-xice0(end));
        if abs(k1) == Inf && k2 ~= Inf
            xtip = xint(i);
            ytip = k2*(xtip - xint(i+1))+yint(i+1);
        elseif k1 ~= Inf && k2 == Inf
            xtip = xint(i+1);
            ytip = k1*(xtip - xint(i))+yint(i);
        elseif k1 == Inf && k2 == Inf
            xtip = nan;ytip = nan;
        else
            xtip = double((k1*xint(i)-yint(i)+yint(i+1)-k2*xint(i+1))/(k1-k2));
            ytip = (abs(k1)<=abs(k2))*(k1*(xtip - xint(i))+yint(i)) + (abs(k1)>abs(k2))*(k2*(xtip - xint(i+1))+yint(i+1));
        end

        aaa = ((xtip-xint(i))^2+(ytip-yint(i))^2)^0.5;
        bbb = ((xtip-xint(i+1))^2+(ytip-yint(i+1))^2)^0.5;
        ccc = ((xint(i)-xint(i+1))^2+(yint(i)-yint(i+1))^2)^0.5;
        if pos(j) == 1 || aaa == 0 || bbb == 0
            the(j) = 180;
        else
            the(j) = acosd((aaa^2+bbb^2-ccc^2)/(2*aaa*bbb));
        end
        if ~isreal(the(j))
%             error('the')
            the(j) = 135;
        end

        alp1 = acosd((aaa^2+ccc^2-bbb^2)/(2*aaa*ccc));
        alp2 = 180-the(j)-alp1;
        alp = max(alp1,alp2);
        posneg = (alp1<=alp2) + (alp1>alp2)*(-1);
        alpha_unknown = alp-(180-the(j))/2;
        alpha_unknown(isnan(alpha_unknown))=0;

        k_DB(j) = tand(alpha_unknown)*posneg;
    end
    icelength = length(iceField.x);
    if meanPhi(j)/pi*180 > 80
        %ice nodes inside contact
        xice_out1 = iceField.x(max(contacts(i)-100,1):contacts(i)-1);
        yice_out1 = iceField.y(max(contacts(i)-100,1):contacts(i)-1);
        xice_out2 = iceField.x(contacts(i+1)+1:min(contacts(i+1)+100,icelength));
        yice_out2 = iceField.y(contacts(i+1)+1:min(contacts(i+1)+100,icelength));
        
        wedge_dist1 = abs((a(j)*xice_out1 + b(j)*yice_out1 + c(j))/sqrt(a(j)^2+b(j)^2));
        wedge_dist2 = abs((a(j)*xice_out2 + b(j)*yice_out2 + c(j))/sqrt(a(j)^2+b(j)^2));
        
        diff1 = diff(wedge_dist1(end:-1:1));
        diff2 = diff(wedge_dist2);
        ind_wedge_dist1 = find(diff1<0,1);
        wedge_dist1 = wedge_dist1(end-ind_wedge_dist1+1);
        ind_wedge_dist2 = find(diff2<0,1);
        wedge_dist2 = wedge_dist2(ind_wedge_dist2);
        
        wedge_end1x = iceField.x(contacts(i)-ind_wedge_dist1);
        wedge_end1y = iceField.y(contacts(i)-ind_wedge_dist1);
        wedge_end2x = iceField.x(contacts(i+1)+ind_wedge_dist2);
        wedge_end2y = iceField.y(contacts(i+1)+ind_wedge_dist2);
                
%         debug
%         scatter(wedge_end1x,wedge_end1y)
%         scatter(wedge_end2x,wedge_end2y)
        if isempty(wedge_dist1) || isempty(wedge_dist2)
            b_wedge(j) = 0;
            c_wedge(j) = 0;
            the(j) = 0;
        else
            b_wedge(j) = mean([wedge_dist1,wedge_dist2]);
            c_wedge(j) = ((xint(i)-xint(i+1))^2+(yint(i)-yint(i+1))^2)^0.5/2;
            a_wedge_minus_c_1 = sqrt( (wedge_end1x-midx(j))^2 + (wedge_end1y-midy(j))^2 - wedge_dist1^2) - c_wedge(j);
            a_wedge_minus_c_2 = sqrt( (wedge_end2x-midx(j))^2 + (wedge_end2y-midy(j))^2 - wedge_dist2^2) - c_wedge(j);
            r_wedge1 = (a_wedge_minus_c_1^2 + wedge_dist1^2) / (2*wedge_dist1);
            r_wedge2 = (a_wedge_minus_c_2^2 + wedge_dist2^2) / (2*wedge_dist2);
            the1 = asind(min(a_wedge_minus_c_1/r_wedge1,1));
            the2 = asind(min(a_wedge_minus_c_2/r_wedge2,1));
            the(j) = 180 - the1 - the2;
        end
        if ~isreal(the(j))
            the(j) = 180;
        end
    end
    %contact area
    area(j) = polyarea(xbound, ybound);
    
    %crushing depth
    if  (contactswl(i) ~= -10000) % -10000 means no ship node inside
        if sign(ship.y(contactswl(i))) ~=sign(ship.y(contactswl(i+1)))
            %                     depth1 = (abs(a(j)*66.1016+c(j))/sqrt(a(j)^2+b(j)^2));
            %                     depth2 = ((xint(i)-xint(i+1))^2+(yint(i)-yint(i+1))^2)^0.5;
            %                     depth(j) = (0.5*depth1*depth2/pi)^0.5;
            depth(j) = isreal(the(j))*(area(j)*2/pi)^0.5;
            LC_DB(j) = 0.33;
            area(j) = area(j)/meanPhi(j);
            ctrx(j) = mean([xint(i),xint(i+1),max(ship.x)]);
            ctry(j) = mean([yint(i),yint(i+1),0]);
        else
            depth(j) = isreal(the(j))*(abs(a(j)*ctrx(j)+b(j)*ctry(j)+c(j))/sqrt(a(j)^2+b(j)^2))/LC_DB(j);
        end
    else
        depth(j) = isreal(the(j))*(abs(a(j)*ctrx(j)+b(j)*ctry(j)+c(j))/sqrt(a(j)^2+b(j)^2))/LC_DB(j);
    end
    
    if meanPhi(j)/pi*180>=80
        ctrx(j) = mean([xint(i),xint(i+1)]);
        ctry(j) = mean([yint(i),yint(i+1)]);
    end
    %
    %             scatter(xice0,yice0);hold on
    %             scatter(xint(i:i+1),yint(i:i+1))
    %             plot(iceField.x,iceField.y,ship.x,ship.y)
    %             xlim([min(xint(i:i+1))-0.1,max(xint(i:i+1))+0.1])
    %             ylim(([min(yint(i:i+1))-0.1,max(yint(i:i+1))+0.1]))
    %             hold off
    %-------------------------------------------------------
    
    midX(j) = r(1) + ctrx(j).*cy3 - ctry(j).*sy3;
    midY(j) = r(2) + ctrx(j).*sy3 + ctry(j).*cy3;
    
    if isnan(midX(j))
        disp('ri le gou')%error
    end
    
%     if sign(ship.y(contactswl_1(i))) ~= sign(ship.y(contactswl_1(i+1)))
%         area(j) = area(j)/meanPhi(j);
%     end
    %             if depth(j)<0
    %                 error('d')
    %             end
    
    if meanPhi(j)/pi*180<80 && depth(j)>0
        %    [Corr_v, Corr_A, Corr_ind] = areacorrec(depth(j),hi(j),sigmac(j),meanPhi(j)/pi*180,iceField.E);
        ind_dep = max(min(round(mean([depth(j)/cosd(the(j)/2),depth(j)])/l_c(j)/0.01)+1,21),1);
        ind_phi = max(round((meanPhi(j)/pi*180-15)/5),1);
        Corr_v = areacorr.V.depth(ind_dep).h.sigc.phi(ind_phi).E;
        Corr_A = areacorr.A.depth(ind_dep).h.sigc.phi(ind_phi).E;
        Corr_ind = areacorr.ind.depth(ind_dep).h.sigc.phi(ind_phi).E;
    else
        Corr_v = 1;
        Corr_A = 1;
        Corr_ind = 1;
    end
    area(j) = area(j)*Corr_A;
    depth(j) = depth(j)*Corr_ind;
%     if meanPhi(j)/pi*180<80 && depth(j)<0
%         disp('debug')
%     end
    j = j + 1;
end

%% Ice bending failure and update ice sheet
%/////////////////////////////////////////////////////////////////////////
%determination of the ice breaking forces. a crushing force normal to thea
%surface of the ship hull is determined first. the crushing force and
%friction forces are then resolved into a vertical and a horizontal
%component


%magnitude of crushing force DETERMINISTINEN JÄÄKENTT? ARVOT

Fcr = crushing(sigmac,area,lin);

%frictional forces
fH = iceField.mu .* Fcr .* abs(vrelt) ./ sqrt(vrelt.^2+vreln1.^2);
fV = iceField.mu .* Fcr .* abs(vreln1) ./ sqrt(vrelt.^2+vreln1.^2);

%crushing forces
FH = Fcr.*smPhi + fV.*cmPhi;%magnitude of FH, direction from normal of contact
FV = Fcr.*cmPhi - fV.*smPhi;%magnitude of FV

%///////////////////////////////////////////////////////////////////
%breaking radius, confined between Rmin and Rmax
%         R = cusp(cra,Cl,Cv,l_c,vreln2,randomNumber);
R = ones(length(vreln2),1)*1;
Redge = R;

% nice = length(iceField.X);%number of ice nodes
distL = 0.05;%discretization length of ice,
dang = distL./R;%incrimental angle (absolute value)
cracks = NaN*ones(ceil(2*pi/min(dang)), 2, length(R));%initialize matrix for X and Y coordinates of crack points
crStEnd = zeros(length(R), 2);
bendCounter = 1;

wedgeAng = nan(length(R),1);

B_ref = [0.005:0.001:0.04,0.042:0.002:0.05,0.055:0.005:0.1];
E = iceField.E;
SF = 0.5*sigmac/1e6.*(E/1e9./hi).^0.5;
SF_H = 0.5*sigmac/1e6.*(E/1e9./hi).^0.25;
rf = 1+1.5*abs(vreln).^0.4;
% rR = (0.6*exp(-10*abs(vreln))+0.4) .* rf.^(-0.1837);
rR = rf.^(-0.1837).*(1-0.1*abs(vreln));

ind_break = [];

x_broken = nan(length(R),2);
for i = 1:length(R)
    
    if meanPhi(i)/pi*180 <= 80 % bow
        nice = length(iceField.X);%number of ice nodes
        %distances from midpoint of contact to ice nodes
        dst = ((midX(i)*ones(1, nice) - iceField.X).^2 + (midY(i)*ones(1, nice) - iceField.Y).^2).^0.5;
        cls = find(dst==min(dst), 1);%index of closest point
        theta = deg2rad(the(i));
        
        if LC_DB(i)>=0.30
            dbdb = DB.DB33;
        elseif LC_DB(i)>=0.24 && LC_DB(i)<0.3
            dbdb = DB.DB27;
        elseif LC_DB(i)>=0.18 && LC_DB(i)<0.24
            dbdb = DB.DB21;
        elseif LC_DB(i)>=0.12 && LC_DB(i)<0.18
            dbdb = DB.DB15;
        elseif LC_DB(i)>=0.06 && LC_DB(i)<0.12
            dbdb = DB.DB9;
        else
            dbdb = DB.DB3;
        end
        
        DBcolumeInd = 1;
        if theta>pi/3 && theta<5/6*pi
            thetaInd = round((theta-pi/3)/(pi/60)+1);
            if meanPhi(i)/pi*180 > 80
                rbreak = Inf;
                Edge1 = 1;
                Edge2 = 1;
                center = 1;
            else
                angleInd = max(round((meanPhi(i)/pi*180-20)/5)+1,1);
                rbreak_ind = zeros(100,1);
                while length(rbreak_ind)==100 || DBcolumeInd==0
                    db = dbdb.DB0;
                    comparedata = db(thetaInd).stress(angleInd,:);
                    comparedataH = db(thetaInd).stressH(angleInd,:);
                    rbreak_ind = find((comparedata*SF(i)+comparedataH*SF_H(i))/rf(i)<sigmaf(i));
                    comparedata(isnan(comparedata))=[];
                    DBcolumeInd = min(length(rbreak_ind)+1,length(comparedata));
                end
                if length(rbreak_ind) == length(comparedata)
                    rbreak = Inf;
                    Edge1 = 1;Edge2 = 1;center = 1;
                else
                    rbreak = B_ref(length(rbreak_ind)+1)*l_c(i);
                    
                    if sign(k_DB(i)) == 1
                        Edge1 = db(thetaInd).edge1(angleInd,DBcolumeInd)*l_c(i)*rR(i);
                        Edge2 = db(thetaInd).edge2(angleInd,DBcolumeInd)*l_c(i)*rR(i);
                    else
                        Edge1 = db(thetaInd).edge2(angleInd,DBcolumeInd)*l_c(i)*rR(i);
                        Edge2 = db(thetaInd).edge1(angleInd,DBcolumeInd)*l_c(i)*rR(i);
                    end
                    center = db(thetaInd).center(angleInd,DBcolumeInd)*l_c(i)*rR(i);
                end
                
            end
            r_min = max( sqrt( (ctrx(i)-xint(2*i-1))^2 + (ctry(i)-yint(2*i-1))^2 ),...
                sqrt( (ctrx(i)-xint(2*i))^2 + (ctry(i)-yint(2*i))^2 ) );
            r_scale = max(r_min/max(Edge1,Edge2),1);
%             if r_scale>1
%                 disp('big cusp')
%             end
            Edge1 = Edge1*r_scale;
            Edge2 = Edge2*r_scale;
            center = center*r_scale;
            %recalculate the angle
            lowInd = cls - find(dst(cls-1:-1:1)>=Edge1, 1);
            highInd = cls + find(dst(cls+1:end)>=Edge2, 1);
            crStEnd(i,:) = [lowInd highInd];
            dxlow = iceField.X(lowInd) - midX(i);
            dylow = iceField.Y(lowInd) - midY(i);
            dxhigh = iceField.X(highInd) - midX(i);
            dyhigh = iceField.Y(highInd) - midY(i);
            lowAng = atan2(dylow,dxlow);
            highAng = atan2(dyhigh,dxhigh);
            %define points along the crack
            if (lowAng>=0&&highAng>=0)||(lowAng<0&&highAng<0)
                N = floor(((highAng>lowAng)*abs(highAng-lowAng) + (highAng<lowAng)*(2*pi-lowAng+highAng))/dang(i));
            elseif lowAng<0&&highAng>=0
                N = floor((abs(lowAng)+highAng)/dang(i));
            elseif lowAng>=0&&highAng<0
                N = floor((2*pi-(abs(highAng)+lowAng))/dang(i));
            end
            
            theta1 = N*dang(i);
            coef = [0 0 1;0.25*theta1^2 0.5*theta1 1;theta1^2 theta1 1]\[Edge1;center;Edge2];
            
            cm = tril(ones(N));%coefficient matrix
            th = lowAng+cm*(dang(i)*ones(N,1));
            rrr = coef(1)*(th-lowAng).^2+coef(2)*(th-lowAng)+coef(3);
            
            cracksx = midX(i)+rrr.*cos(th);
            cracksy = midY(i)+rrr.*sin(th);
        else
            db = dbdb.DB0;
            thetaInd = 1;
            if meanPhi(i)/pi*180 > 80
                rbreak = Inf;
                rrr = 1;
            else
                angleInd = max(round((meanPhi(i)/pi*180-20)/5)+1,1);
                rbreak_ind = find((db(thetaInd).stress(angleInd,:)*SF(i)+ ...
                    db(thetaInd).stressH(angleInd,:)*SF_H(i))/rf(i)<sigmaf(i));
                
                if length(rbreak_ind) == 51
                    rbreak = Inf;
                else
                    rbreak = B_ref(length(rbreak_ind)+1)*l_c(i);
                end
                DBcolumeInd = min(length(rbreak_ind)+1,51);
                
                rrr = db(thetaInd).center(angleInd,DBcolumeInd)*l_c(i)*rR(i);
            end
            r_min = max( sqrt( (ctrx(i)-xint(2*i-1))^2 + (ctry(i)-yint(2*i-1))^2 ),...
                sqrt( (ctrx(i)-xint(2*i))^2 + (ctry(i)-yint(2*i))^2 ) );
            r_scale = max(r_min/max(rrr),1);
            rrr = rrr*r_scale;
            %recalculate the angle
            lowInd = cls - find(dst(cls-1:-1:1)>=rrr, 1);
            highInd = cls + find(dst(cls+1:end)>=rrr, 1);
            crStEnd(i,:) = [lowInd highInd];
            dxlow = iceField.X(lowInd) - midX(i);
            dylow = iceField.Y(lowInd) - midY(i);
            dxhigh = iceField.X(highInd) - midX(i);
            dyhigh = iceField.Y(highInd) - midY(i);
            lowAng = atan2(dylow,dxlow);
            highAng = atan2(dyhigh,dxhigh);
            %define points along the crack
            if (lowAng>=0&&highAng>=0)||(lowAng<0&&highAng<0)
                N = floor(((highAng>lowAng)*abs(highAng-lowAng) + (highAng<lowAng)*(2*pi-lowAng+highAng))/dang(i));
            elseif lowAng<0&&highAng>=0
                N = floor((abs(lowAng)+highAng)/dang(i));
            elseif lowAng>=0&&highAng<0
                N = floor((2*pi-(abs(highAng)+lowAng))/dang(i));
            end
            
            cm = tril(ones(N));%coefficient matrix
            th = lowAng+cm*(dang(i)*ones(N,1));
            
            cracksx = midX(i)+rrr*cos(th);
            cracksy = midY(i)+rrr*sin(th);
            
        end
        
        
        if lowInd == []
            disp('ri')
        end
        
        
        cracksy(cracksx<min(iceField.X)) = [];
        cracksx(cracksx<min(iceField.X)) = [];
        %limiting force for bending failure
        %             Pf = 1*bearing(dyna,N*dang(i),sigmaf,hi(i),vreln2(i),Cf);

        %check for bending
        if (depth(i)>=rbreak && ~( (r_scale>1 || cLength(i)>0.5*l_c(i)) && midx(i)<ship.mids)) || (cLength(i) > 0.5*l_c(i) && midx(i)>ship.mide)
            
            if max(rrr) > l_c(i)
                disp('big cusp')
            end
            ind_break(end+1) = i;
            %                disp([num2str(rbreak),'    ',num2str(area(i)),'    ',num2str(Fcr(i))])
            %                 ratio = Pf/FV(i);
            %                 FV(i) = Pf;
            %                 Fcr(i) = ratio*Fcr(i);
            %                 fV(i) = (iceField.mu*Fcr(i)*abs(vreln1(i)))/norm([vrelt(i) vreln1(i)]);
            %                 FH(i) = Fcr(i)*smPhi(i) + fV(i)*cmPhi(i);
            %                 fH(i) = iceField.mu.*Fcr(i).*abs(vrelt(i))./norm([vrelt(i) vreln1(i)]);
            bendCounter = bendCounter+1;
            lowInd = cls - find(dst(cls-1:-1:1)>=rrr(1), 1);
            highInd = cls + find(dst(cls+1:end)>=rrr(end), 1);
            
            %                 plot(iceField.X(1:lowInd),iceField.Y(1:lowInd))
            %                 hold on
            %                 plot(iceField.X(highInd:end),iceField.Y(highInd:end))
            %                 scatter(cracksx,cracksy)
            %                 if isempty(cracksx)
            %                     disp('error')
            %                 end
            %                 xlim([min(cracksx)-1,max(cracksx)+1])
            %                 ylim([min(cracksy)-1,max(cracksy)+1])
            %                 hold off
            
            iceField.X = horzcat(iceField.X(1:lowInd), cracksx', iceField.X(highInd:end));
            iceField.Y = horzcat(iceField.Y(1:lowInd), cracksy', iceField.Y(highInd:end));
            nice = length(iceField.X);
            wedgeAng(i) = N*dang(i);
            if length(rrr)==1
                R(i) = rrr;
                Redge(i) = rrr;
            else
                R(i) = center;
                Redge(i) = (Edge1 + Edge2)*0.5;
            end
            icepiece = iceField.icepiece;
            if the(i) == 180
                icepiecenum = 3;
                thetaw = deg2rad(80);
                icepiece.loc(icepiece.num+[1:3],:) = [midx(i),0.001;midx(i),-0.001;midx(i),0];
            elseif the(i) == 0
                icepiecenum = 0;
            else
                icepiecenum = 1;
                thetaw = theta;
                icepiece.loc(icepiece.num+1,:) = [midx(i),midy(i)];
            end
            for iii = 1:icepiecenum
                icepiece.num = icepiece.num + 1;
                icepiece.clength(icepiece.num) = cLength(i); % contact length
                icepiece.RelDis(icepiece.num) = 0; 
                if atan(hi(i)/R(i))+meanPhi(i)>90/180*pi % not rotatable
                    icepiece.IndRot(icepiece.num) = 0; 
                else                                        % rotatable
                    icepiece.IndRot(icepiece.num) = 1; 
                end
%                icepiece.loc(icepiece.num,:) = [midx(i),midy(i)];
%                icepiece.theta(icepiece.num) = 0;
                icepiece.L(icepiece.num) = R(i);

                icepiece.Lprime(icepiece.num) = Redge(i)*cos(thetaw/2);
                icepiece.thetaw(icepiece.num) = thetaw;
                k1 = 2*tan(thetaw/2);
                k2 = 2*icepiece.Lprime(icepiece.num)*tan(thetaw/2)/(icepiece.L(icepiece.num)-icepiece.Lprime(icepiece.num));

%                 icepiece.Lprime(icepiece.num) = Redge(i)*sin(theta/2);
%                 icepiece.thetaw(icepiece.num) = theta;
%                 k1 = 2*tan(theta/2);
%                 k2 = 2*tan(icepiece.Lprime*tan(theta/2))/(icepiece.L-icepiece.Lprime);

                icepiece.S(icepiece.num) = icepiece.Lprime(icepiece.num)*k1*icepiece.Lprime(icepiece.num)/2 + ...
                    (icepiece.L(icepiece.num)-icepiece.Lprime(icepiece.num))*k2*(icepiece.L(icepiece.num)-icepiece.Lprime(icepiece.num))/2;
                icepiece.Sc(icepiece.num) = k1*(1/2*icepiece.L(icepiece.num)*icepiece.Lprime(icepiece.num)^2-1/3*icepiece.Lprime(icepiece.num)^3) +...
                    1/3*k2*(icepiece.L(icepiece.num)-icepiece.Lprime(icepiece.num))^3;
                if icepiece.Sc(icepiece.num)<0
                    disp('Sc less than 0')
                end
                if isnan(icepiece.Sc(icepiece.num)) ||isnan(icepiece.S(icepiece.num))
                    error('Sc is nan')
                end
                icepiece.SrhoBsq(icepiece.num) = k1*(1/2*icepiece.L(icepiece.num)^2*icepiece.Lprime(icepiece.num)^2-2/3*icepiece.L(icepiece.num)*icepiece.Lprime(icepiece.num)^3+...
                    1/4*icepiece.Lprime(icepiece.num)^4) + 1/4*k2*(icepiece.L(icepiece.num)-icepiece.Lprime(icepiece.num))^4;
                icepiece.stage(icepiece.num) = 0;
                icepiece.vn(icepiece.num) = abs(vreln(i));
                icepiece.Frotate = zeros(4,1);
                icepiece.x(icepiece.num,:) = [0,0,0];
                icepiece.y(icepiece.num,:) = [0,0,0];
                icepiece.h(icepiece.num) = hi(i);
                icepiece.t_initial(icepiece.num) = 0;
                icepiece.F_initial(icepiece.num) = 0;
                icepiece.Fn(icepiece.num) = 0;
                icepiece.theta(icepiece.num,:) = [0,0,0];
    %             icepiece.phi = meanPhi(i);
    %             icepiece.mu_n(icepiece.num) = iceField.mu*vreln(i)*tan(meanPhi(i))/sqrt(vrelt(i)^2+(vreln(i)*tan(meanPhi(i)))^2);
                icepiece.s(icepiece.num) = 0;
                icepiece.pos(icepiece.num) = 0;
                icepiece.m(icepiece.num,:) = [hi(i)*icepiece.S(icepiece.num)*iceField.rhoi;...
                    hi(i)*icepiece.S(icepiece.num)*iceField.rhoi;...
                    icepiece.SrhoBsq(icepiece.num)*hi(i)*iceField.rhoi-hi(i)*icepiece.S(icepiece.num)*iceField.rhoi*(icepiece.Sc(icepiece.num)/icepiece.S(icepiece.num))^2];
            end
            if icepiecenum == 3
                icepiece.pos(icepiece.num) = 1;
            end
            iceField.icepiece = icepiece;
        else
            cracks(:,:,bendCounter) = NaN;
            R(i) = NaN;
            Redge(i) = NaN;
            crStEnd(i) = NaN;
        end
    else % midship
        %        xintsort = sort([xint(2*i-1),xint(2*i)]);
        %        if (xintsort(1)>=ship.mids && xintsort(2)<=ship.mide)% || (xintsort(2)>=ship.mids && xintsort(2)<=ship.mide)
        %             [sep_plus, sep_minus] = wedgeseparation(iceField.x,iceField.y,[ship.mids,ship.mide]);
        %
        %             if yint(2*i-1)<0
        %                 sorted = sortrows(sep_minus,1);
        %             else
        %                 sorted = sortrows(sep_plus,1);
        %             end
        %             pos1 = find((sorted(:,1)-xintsort(1))<=0);
        %             pos1 = pos1(end);
        %             pos2 = find((sorted(:,1)-xintsort(2))>=0);
        %             pos2 = pos2(1);
        %             y_sep1 = sorted(pos1,2);
        %             y_sep2 = sorted(pos2,2);
        
        %             b_wedge(i) = abs(abs(mean([y_sep1,y_sep2]))-abs(yint(2*i-1)))/l_c(i);
        %             c_wedge(i) = diff(xintsort)/l_c(i)/2;
        %h and the were calculated previously
        theta = deg2rad(the(i));
        scale_midship = (E/hi(i))^0.25 / (50e6/0.03)^0.25; % needs correction
        %  fprintf("angle=%f,b=%f,c=%f,h=%f,mu=%f\n",the(i),b_wedge(i),c_wedge(i),hi(i),abs(FV(i)/FH(i)))
        %         if b_wedge(i)>10
        %             pause(0.1)
        %         end
        sigma_midship = 0.5*sigmac(i)*NN1([the(i),b_wedge(i)/l_c(i),c_wedge(i)/l_c(i),abs(FV(i)/FH(i)),abs(fH(i)/FH(i))]);
        R_midship = rR(i)*l_c(i)*NN2([the(i),b_wedge(i)/l_c(i),c_wedge(i)/l_c(i),abs(FV(i)/FH(i)),abs(fH(i)/FH(i))]);
        %  plot(iceField.x,iceField.y,ship.x,ship.y)
        sigma_midship = sigma_midship*scale_midship;
%         if abs(FV(i)/FH(i))>0.15 && c_wedge(i)/l_c(i)>0.05
%             pause(0.1)
%         end
         if sigma_midship/rf(i)>iceField.meansf && c_wedge(i)/l_c(i)<1
%          if 1<0
            if max(R_midship) > l_c(i)
                disp('big cusp')
            end
            
            ind_break(end+1) = i;
            x_broken(i,:) = [midx(i)-R_midship, midx(i)+R_midship];
            nice = length(iceField.X);
            dst = ((midX(i)*ones(1, nice) - iceField.X).^2 + (midY(i)*ones(1, nice) - iceField.Y).^2).^0.5;
            cls = find(dst==min(dst), 1);%index of closest point
            lowInd = cls - find(dst(cls-1:-1:1)>=R_midship, 1);
            highInd = cls + find(dst(cls+1:end)>=R_midship, 1);
            dxlow = iceField.X(lowInd) - midX(i);
            dylow = iceField.Y(lowInd) - midY(i);
            dxhigh = iceField.X(highInd) - midX(i);
            dyhigh = iceField.Y(highInd) - midY(i);
            lowAng = atan2(dylow,dxlow);
            highAng = atan2(dyhigh,dxhigh);
            %define points along the crack
            if (lowAng>=0&&highAng>=0)||(lowAng<0&&highAng<0)
                N = floor(((highAng>lowAng)*abs(highAng-lowAng) + (highAng<lowAng)*(2*pi-lowAng+highAng))/dang(i));
            elseif lowAng<0&&highAng>=0
                N = floor((abs(lowAng)+highAng)/dang(i));
            elseif lowAng>=0&&highAng<0
                N = floor((2*pi-(abs(highAng)+lowAng))/dang(i));
            end
            
%             cm = tril(ones(N));%coefficient matrix
%             th = lowAng+cm*(dang(i)*ones(N,1));
%            
            theta1 = N*dang(i);
            coef = [0 0 1;0.25*theta1^2 0.5*theta1 1;theta1^2 theta1 1]\[R_midship;mean([R_midship,R_midship*cosd(the(i)/2)]);R_midship];
            
            cm = tril(ones(N));%coefficient matrix
            th = lowAng+cm*(dang(i)*ones(N,1));
            rrr = coef(1)*(th-lowAng).^2+coef(2)*(th-lowAng)+coef(3);
            
            cracksx = midX(i)+rrr.*cos(th);
            cracksy = midY(i)+rrr.*sin(th);
            
            lowInd = cls - find(dst(cls-1:-1:1)>=R_midship, 1);
            highInd = cls + find(dst(cls+1:end)>=R_midship, 1);
            
            iceField.X = horzcat(iceField.X(1:lowInd), cracksx', iceField.X(highInd:end));
            iceField.Y = horzcat(iceField.Y(1:lowInd), cracksy', iceField.Y(highInd:end));
           
            % Turning
            icepiece = iceField.icepiece;
            icepiecenum = 1;
            thetaw = theta;
            icepiece.loc(icepiece.num+1,:) = [midx(i),midy(i)];

            for iii = 1:icepiecenum
                icepiece.num = icepiece.num + 1;
                icepiece.clength(icepiece.num) = cLength(i);
%                 icepiece.loc(icepiece.num,:) = [midx(i),midy(i)];
%             icepiece.theta(icepiece.num) = 0;
                icepiece.L(icepiece.num) = mean([R_midship,R_midship*cosd(the(i)/2)]);
                icepiece.RelDis(icepiece.num) = 0; 
                if atan(hi(i)/icepiece.L(icepiece.num))+meanPhi(i)>90/180*pi % not rotatable
                    icepiece.IndRot(icepiece.num) = 0; 
                else                                        % rotatable
                    icepiece.IndRot(icepiece.num) = 1; 
                end
                icepiece.Lprime(icepiece.num) = R_midship*cos(thetaw/2);
                icepiece.thetaw(icepiece.num) = thetaw;
                k1 = 2*tan(thetaw/2);
                k2 = 2*icepiece.Lprime(icepiece.num)*tan(thetaw/2)/(icepiece.L(icepiece.num)-icepiece.Lprime(icepiece.num));

%                 icepiece.Lprime(icepiece.num) = Redge(i)*sin(theta/2);
%                 icepiece.thetaw(icepiece.num) = theta;
%                 k1 = 2*tan(theta/2);
%                 k2 = 2*tan(icepiece.Lprime*tan(theta/2))/(icepiece.L-icepiece.Lprime);

                icepiece.S(icepiece.num) = icepiece.Lprime(icepiece.num)*k1*icepiece.Lprime(icepiece.num)/2 + ...
                    (icepiece.L(icepiece.num)-icepiece.Lprime(icepiece.num))*k2*(icepiece.L(icepiece.num)-icepiece.Lprime(icepiece.num))/2;
                icepiece.Sc(icepiece.num) = k1*(1/2*icepiece.L(icepiece.num)*icepiece.Lprime(icepiece.num)^2-1/3*icepiece.Lprime(icepiece.num)^3) +...
                    1/3*k2*(icepiece.L(icepiece.num)-icepiece.Lprime(icepiece.num))^3;
                if icepiece.Sc(icepiece.num)<0
                    error('Sc less than 0')
                end
                if isnan(icepiece.Sc(icepiece.num)) || isnan(icepiece.S(icepiece.num))
                    error('Sc is nan')
                end
                icepiece.SrhoBsq(icepiece.num) = k1*(1/2*icepiece.L(icepiece.num)^2*icepiece.Lprime(icepiece.num)^2-2/3*icepiece.L(icepiece.num)*icepiece.Lprime(icepiece.num)^3+...
                    1/4*icepiece.Lprime(icepiece.num)^4) + 1/4*k2*(icepiece.L(icepiece.num)-icepiece.Lprime(icepiece.num))^4;
                icepiece.stage(icepiece.num) = 0;
                icepiece.vn(icepiece.num) = abs(vreln(i));
                icepiece.Frotate = zeros(4,1);
                icepiece.x(icepiece.num,:) = [0,0,0];
                icepiece.y(icepiece.num,:) = [0,0,0];
                icepiece.h(icepiece.num) = hi(i);
                icepiece.t_initial(icepiece.num) = 0;
                icepiece.F_initial(icepiece.num) = 0;
                icepiece.Fn(icepiece.num) = 0;
                icepiece.theta(icepiece.num,:) = [0,0,0];
    %             icepiece.phi = meanPhi(i);
    %             icepiece.mu_n(icepiece.num) = iceField.mu*vreln(i)*tan(meanPhi(i))/sqrt(vrelt(i)^2+(vreln(i)*tan(meanPhi(i)))^2);
                icepiece.s(icepiece.num) = 0;
                icepiece.pos(icepiece.num) = 0;
                icepiece.m(icepiece.num,:) = [hi(i)*icepiece.S(icepiece.num)*iceField.rhoi;...
                    hi(i)*icepiece.S(icepiece.num)*iceField.rhoi;...
                    icepiece.SrhoBsq(icepiece.num)*hi(i)*iceField.rhoi-hi(i)*icepiece.S(icepiece.num)*iceField.rhoi*(icepiece.Sc(icepiece.num)/icepiece.S(icepiece.num))^2];
            end
            iceField.icepiece = icepiece;
%             iceField.x = 1/(cos(r(6))^2+sin(r(6))^2)*(iceField.X*cos(r(6))+iceField.Y*sin(r(6)))-(r(1)*cos(r(6))+r(2)*sin(r(6)));
%             iceField.y = 1/(cos(r(6))^2+sin(r(6))^2)*(-iceField.X*sin(r(6))+iceField.Y*cos(r(6)))+(r(1)*sin(r(6))-r(2)*cos(r(6)));
        else
%             cracks(:,:,bendCounter) = NaN;
            R(i) = NaN;
            Redge(i) = NaN;
            crStEnd(i) = NaN;
        end
    end
end
x_broken(isnan(x_broken(:,1)),:) = [];
%% Crushing at midship area
iceField.x = (iceField.X*cos(r(4))+iceField.Y*sin(r(4)))-(r(1)*cos(r(4))+r(2)*sin(r(4)));
iceField.y = (-iceField.X*sin(r(4))+iceField.Y*cos(r(4)))+(r(1)*sin(r(4))-r(2)*cos(r(4)));
for i = 1:length(xint)/2
    unbroken = 1;
    if abs(xint(2*i-1)-xint(2*i))>0.02
        [xintsort,ind] = sort([xint(2*i-1),xint(2*i)]);
        yintsort = [yint(2*i-1),yint(2*i)];
        yintsort = yintsort(ind);
        for j = 1:length(x_broken(:,1))
            if ~(min(x_broken(j,:))>max(xintsort) || max(x_broken(j,:))<min(xintsort))
                unbroken = 0*unbroken;
            end
        end
        if (xintsort(2)>=ship.mids && xintsort(2)<=ship.mide) && isempty(find(ind_break == i, 1)) && unbroken == 1% || (xintsort(2)>=ship.mids && xintsort(2)<=ship.mide)
            ind_hull = ( ship.x > xintsort(1) & ship.x < xintsort(2) & sign(ship.y) == sign(yint(2*i)) );
            crackx_hull_unsorted = ship.x(ind_hull);
            cracky_hull_unsorted = ship.y(ind_hull);
            [crackx_hull,ind] = sort(crackx_hull_unsorted);
            cracky_hull = cracky_hull_unsorted(ind);
            
            cracksx_interp = xintsort(1):0.1:xintsort(2);
            cracksx_interp([1,end]) = [];
            
            if ~isempty(crackx_hull_unsorted)
                cracksy_interp = interp1([xintsort(1),crackx_hull',xintsort(2)],[yintsort(1),cracky_hull',yintsort(2)],cracksx_interp);
            else
                cracksx_interp = [];
                cracksy_interp = [];
            end
            
            if xintsort(1) == xint(2*i-1)
                cracksxmid = [xint(2*i-1),cracksx_interp,xint(2*i)];
                cracksymid = [yint(2*i-1),cracksy_interp,yint(2*i)];
            else
                cracksxmid = [xint(2*i-1),cracksx_interp(end:-1:1),xint(2*i)];
                cracksymid = [yint(2*i-1),cracksy_interp(end:-1:1),yint(2*i)];
            end
            
            while sum(abs(diff(cracksxmid))<0.01)>0
                diffcrack = diff(cracksxmid);
                ind_del = find(diffcrack<0.01,1)+1;
                if ind_del < length(cracksxmid)
                    cracksxmid(ind_del) = [];
                    cracksymid(ind_del) = [];
                else
                    break
                end
            end

            cracksXmid = r(1) + cracksxmid*cy3 - cracksymid*sy3;
            cracksYmid = r(2) + cracksxmid*sy3 + cracksymid*cy3;
            if yint(2*i-1)<0
                lowIndmidall = find(iceField.x(1:end-4)<xint(2*i-1) & iceField.y(1:end-4)<0 & iceField.x(1:end-4)>ship.min & abs(iceField.y(1:end-4))<ship.B/2+2*l_c(i)); %...
                  %   & iceField.x(1:end-4)<ship.max & abs(iceField.y(1:end-4))<ship.B/2+2*l_c(i));
                lowIndmid = lowIndmidall(end);
                ind_minus = 1;
                while abs(iceField.x(lowIndmid)-xint(2*i-1))<0.01
                    lowIndmid = lowIndmidall(end-ind_minus);
                    ind_minus = ind_minus+1;
                end
                highIndmid = find(iceField.x(1:end-4)>xint(2*i) & iceField.y(1:end-4)<=ship.B/2 & iceField.x(1:end-4)>ship.min & abs(iceField.y(1:end-4))<ship.B/2+2*l_c(i),1); % ...
                  %  & iceField.x(1:end-4)>ship.min & iceField.x(1:end-4)<ship.max & abs(iceField.y(1:end-4))<ship.B/2+2*l_c(i),1);
            else
                lowIndmidall = find(iceField.x(1:end-4)>xint(2*i-1) & iceField.y(1:end-4)>=ship.B/2 & iceField.x(1:end-4)>ship.min & abs(iceField.y(1:end-4))<ship.B/2+2*l_c(i)); % ...
                  %  & iceField.x(1:end-4)>ship.min & iceField.x(1:end-4)<ship.max & abs(iceField.y(1:end-4))<ship.B/2+2*l_c(i));
                lowIndmid = lowIndmidall(end);
                highIndmidall = find(iceField.x(1:end-4)<xint(2*i) & iceField.y(1:end-4)>0 & iceField.x(1:end-4)>ship.min & abs(iceField.y(1:end-4))<ship.B/2+2*l_c(i)); % ...
                  %  & iceField.x(1:end-4)>ship.min & iceField.x(1:end-4)<ship.max & abs(iceField.y(1:end-4))<ship.B/2+2*l_c(i));
                highIndmid = highIndmidall(1);
                ind_minus = 1;
                while abs(iceField.x(highIndmid)-xint(2*i))<0.01
                    highIndmid = highIndmidall(1+ind_minus);
                    ind_minus = ind_minus+1;
                end
            end
            if isempty(lowIndmid) || isempty(lowIndmid)
                error('error in crushing')
            end
            
%             if length(cracksXmid) >100
%                 pause(0.1)
%             end
            ice_debug = iceField;
            iceField.X = horzcat(iceField.X(1:lowIndmid), cracksXmid, iceField.X(highIndmid:end));
            iceField.Y = horzcat(iceField.Y(1:lowIndmid), cracksYmid, iceField.Y(highIndmid:end));
            iceField.x = horzcat(iceField.x(1:lowIndmid), cracksxmid, iceField.x(highIndmid:end));
            iceField.y = horzcat(iceField.y(1:lowIndmid), cracksymid, iceField.y(highIndmid:end));
            
            if sum(abs(diff(iceField.x(1:end-4)))>0.5)>0
                plot(iceField.x,iceField.y,ship.x,ship.y)
                hold on
                scatter(xint,yint)
                hold off
            end
        end
    end
end


% iceField.X = round(iceField.X,3);
% iceField.Y = round(iceField.Y,3);
% iceField.x = round(iceField.x,3);
% iceField.y = round(iceField.y,3);

ind = (diff(iceField.X) == 0) & (diff(iceField.Y) == 0);
if sum(ind)>0
    iceField.X(ind) = [];
    iceField.Y(ind) = [];
    iceField.x(ind) = [];
    iceField.y(ind) = [];
end

%% /////////////////////////////////////////////////////////////////////////
%/////////////////////////////////////////////////////////////////////////
%forces in the x- and y direction
%         if jjj ~= 1
xi = repmat([1 0], length(fH), 1);
yi = repmat([0 1], length(fH), 1);
Fx = dot(tangent,xi,2).*fH + dot(normal,xi,2).*FH;
Fy = dot(tangent,yi,2).*fH + dot(normal,yi,2).*FH;
% Fz = FV;
% Mx = Fz.*midy(:)-Fy.*ship.zarm;
Mx = -Fy.*ship.zarm;
% My = -Fz.*midx(:)+Fx.*ship.zarm;
Mz = -Fx .* midy(:) + Fy .* midx(:);
% compute total x and y direction forces and moment
% iceForce = [sum(Fx(isfinite(Fx))); sum(Fy(isfinite(Fy))); sum(Fz(isfinite(Fz))); sum(Mx(isfinite(Mx))); sum(My(isfinite(My))); sum(Mz(isfinite(Mz)))];
iceForce = [sum(Fx(isfinite(Fx))); sum(Fy(isfinite(Fy))); sum(Mx(isfinite(Mx))); sum(Mz(isfinite(Mz)))];
k_DB(isempty(k_DB)) = nan;
nLoad = [midx(:) midy(:) cLength(:) Fcr(:) FH(:) FV(:) fH(:) fV(:) LC_DB(:) sigmac(:) sigmaf(:) hi(:) R(:) Redge(:) the(:) depth(:)];
end

% OLD VERSION
% for i = 1:length(xint)/2
%     unbroken = 1;
%     if abs(xint(2*i-1)-xint(2*i))>0.00001
%         xintsort = sort([xint(2*i-1),xint(2*i)]);
%         for j = 1:length(x_broken(:,1))
%             if ~(min(x_broken(j,:))>max(xintsort) || max(x_broken(j,:))<min(xintsort))
%                 unbroken = 0*unbroken;
%             end
%         end
%         if (xintsort(1)>=ship.mids && xintsort(2)<=ship.mide) && isempty(find(ind_break == i, 1)) && unbroken == 1% || (xintsort(2)>=ship.mids && xintsort(2)<=ship.mide)
% 
%             n = max(ceil(abs(xint(2*i)-xint(2*i-1))/0.1),1);%number of nodes of new crack
%             if diff(xintsort)>0.1
%                 cracksxmid = xint(2*i-1):(xint(2*i)-xint(2*i-1))/n:xint(2*i);
%                 cracksymid = yint(2*i-1)*ones(1,n+1);
%             else
%                 cracksxmid = [xint(2*i-1),xint(2*i)];
%                 cracksymid = yint(2*i-1)*ones(1,2);
%             end
%             cracksXmid = r(1) + cracksxmid*cy3 - cracksymid*sy3;
%             cracksYmid = r(2) + cracksxmid*sy3 + cracksymid*cy3;
%             if yint(2*i-1)<0
%                 lowIndmidall = find(iceField.x(1:end-4)<xint(2*i-1)&iceField.y(1:end-4)<yint(2*i-1)...
%                     &iceField.x(1:end-4)>ship.min&iceField.x(1:end-4)<ship.max&abs(iceField.y(1:end-4))<ship.B/2+2*l_c(i));
%                 lowIndmid = lowIndmidall(end);
%                 if isempty(lowIndmidall)
% %                     closest = min(iceField.x(iceField.y<yint(2*i-1)) - xint(2*i-1))+xint(2*i-1);
% %                     lowIndmidall = find(iceField.x == closest);
% %                     lowIndmid = find(iceField.x == min(iceField.x(iceField.y<-ship.B/2*2)));
%                     lowIndmidall = find(iceField.x(1:end-4)<xint(2*i-1)&iceField.y(1:end-4)<0 ...
%                         &iceField.x(1:end-4)>ship.min&iceField.x(1:end-4)<ship.max&abs(iceField.y(1:end-4))<ship.B/2+2*l_c(i));
%                     lowIndmid = lowIndmidall(end);
%                 end
%                 highIndmid = find(iceField.x(1:end-4)>xint(2*i)&iceField.y(1:end-4)<yint(2*i-1)...
%                     &iceField.x(1:end-4)>ship.min&iceField.x(1:end-4)<ship.max&abs(iceField.y(1:end-4))<ship.B/2+2*l_c(i),1);
%             else
%                 lowIndmidall = find(iceField.x(1:end-4)>xint(2*i-1)&iceField.y(1:end-4)>yint(2*i-1)...
%                     &iceField.x(1:end-4)>ship.min&iceField.x(1:end-4)<ship.max&abs(iceField.y(1:end-4))<ship.B/2+2*l_c(i));
%                 lowIndmid = lowIndmidall(end);
%                 highIndmid = find(iceField.x(1:end-4)<xint(2*i)&iceField.y(1:end-4)>yint(2*i-1)...
%                     &iceField.x(1:end-4)>ship.min&iceField.x(1:end-4)<ship.max&abs(iceField.y(1:end-4))<ship.B/2+2*l_c(i),1);
%                 if isempty(highIndmid)
% %                     closest = min(iceField.x(iceField.y>yint(2*i-1)) - xint(2*i-1))+xint(2*i-1);
% %                     highIndmid1 = find(iceField.x == closest);
% %                     highIndmid = highIndmid1(1);
% %                     highIndmid = find(iceField.x == min(iceField.x(iceField.y>ship.B/2*2)));
%                     highIndmid = find(iceField.x(1:end-4)<xint(2*i)&iceField.y(1:end-4)>0 ...
%                         &iceField.x(1:end-4)>ship.min&iceField.x(1:end-4)<ship.max&abs(iceField.y(1:end-4))<ship.B/2+2*l_c(i),1);
%                 end
%             end
%             if isempty(lowIndmid) || isempty(lowIndmid)
%                 error('error in crushing')
%             end
%             
%             if length(cracksXmid) >100
%                 pause(0.1)
%             end
%             
%             iceField.X = horzcat(iceField.X(1:lowIndmid), cracksXmid, iceField.X(highIndmid:end));
%             iceField.Y = horzcat(iceField.Y(1:lowIndmid), cracksYmid, iceField.Y(highIndmid:end));
%             iceField.x = horzcat(iceField.x(1:lowIndmid), cracksxmid, iceField.x(highIndmid:end));
%             iceField.y = horzcat(iceField.y(1:lowIndmid), cracksymid, iceField.y(highIndmid:end));
%         elseif (xintsort(1)<=ship.mids && xintsort(2)>ship.mids) && isempty(find(ind_break == i, 1)) && unbroken == 1
%             n1 = max(ceil(abs(xintsort(2)-ship.mids)/0.1),1);%number of nodes of new crack
%             n2 = max(ceil(abs(xintsort(1)-ship.mids)/0.1),1);%number of nodes of new crack
%             if xintsort(2) == xint(2*i-1)
%                 if diff(xintsort)>0.1
%                     cracksxmid = [xint(2*i-1):(ship.mids-xint(2*i-1))/n1:ship.mids,ship.mids:(xint(2*i)-ship.mids)/n2:xint(2*i)];
%                     cracksymid = [yint(2*i-1)*ones(1,n1+1),yint(2*i-1):(yint(2*i)-yint(2*i-1))/n2:yint(2*i)];
%                 else
%                     cracksxmid = [xint(2*i-1),ship.mids,xint(2*i)];
%                     cracksymid = [yint(2*i-1)*ones(1,2),yint(2*i)];
%                 end
%             else
%                 if diff(xintsort)>0.1
%                     cracksxmid = [xint(2*i-1):(ship.mids-xint(2*i-1))/n2:ship.mids,ship.mids:(xint(2*i)-ship.mids)/n1:xint(2*i)];
%                     cracksymid = [yint(2*i-1):(yint(2*i)-yint(2*i-1))/n2:yint(2*i),yint(2*i)*ones(1,n1+1)];
%                 else
%                     cracksxmid = [xint(2*i-1),ship.mids,xint(2*i)];
%                     cracksymid = [yint(2*i-1),yint(2*i)*ones(1,2)];
%                 end
%             end
%             cracksXmid = r(1) + cracksxmid*cy3 - cracksymid*sy3;
%             cracksYmid = r(2) + cracksxmid*sy3 + cracksymid*cy3;
% %             if xintsort(2) == xint(2*i-1)
%             if yint(2*i-1)<0
%                 lowIndmidall = find(iceField.x(1:end-4)<xint(2*i-1)&iceField.y(1:end-4)<0 ...
%                     &iceField.x(1:end-4)>ship.min&iceField.x(1:end-4)<ship.max&abs(iceField.y(1:end-4))<ship.B/2+2*l_c(i));
%                 lowIndmid = lowIndmidall(end);
% %                     if isempty(lowIndmidall)
% %     %                     closest = min(iceField.x(iceField.y<yint(2*i-1)) - xint(2*i-1))+xint(2*i-1);
% %     %                     lowIndmidall = find(iceField.x == closest);
% %     %                     lowIndmid = find(iceField.x == min(iceField.x(iceField.y<-ship.B/2*2)));
% %                         lowIndmidall = find(iceField.x(1:end-4)<xint(2*i-1)&iceField.y(1:end-4)<yint(2*i-1)+0.1...
% %                             &iceField.x(1:end-4)>ship.min&iceField.x(1:end-4)<ship.max&abs(iceField.y(1:end-4))<ship.B/2+2*l_c(i));
% %                         lowIndmid = lowIndmidall(end);
% %                     end
%                 highIndmid = find(iceField.x(1:end-4)>xint(2*i) & iceField.y(1:end-4)<0 ...
%                     &iceField.x(1:end-4)>ship.min & iceField.x(1:end-4)<ship.max & abs(iceField.y(1:end-4))<ship.B/2+2*l_c(i),1);
%             else
%                 lowIndmidall = find(iceField.x(1:end-4)>xint(2*i-1)&iceField.y(1:end-4)>0 ...
%                     &iceField.x(1:end-4)>ship.min&iceField.x(1:end-4)<ship.max&abs(iceField.y(1:end-4))<ship.B/2+2*l_c(i));
%                 lowIndmid = lowIndmidall(end);
%                 highIndmid = find(iceField.x(1:end-4)<xint(2*i)&iceField.y(1:end-4)>0 ...
%                     &iceField.x(1:end-4)>ship.min&iceField.x(1:end-4)<ship.max&abs(iceField.y(1:end-4))<ship.B/2+2*l_c(i),1);
% %                     if isempty(highIndmid)
% %     %                     closest = min(iceField.x(iceField.y>yint(2*i-1)) - xint(2*i-1))+xint(2*i-1);
% %     %                     highIndmid1 = find(iceField.x == closest);
% %     %                     highIndmid = highIndmid1(1);
% %     %                     highIndmid = find(iceField.x == min(iceField.x(iceField.y>ship.B/2*2)));
% %                         highIndmid = find(iceField.x(1:end-4)<xint(2*i) & iceField.y(1:end-4)>yint(2*i)-0.1 ...
% %                             &iceField.x(1:end-4)>ship.min&iceField.x(1:end-4)<ship.max&abs(iceField.y(1:end-4))<ship.B/2+2*l_c(i),1);
% %                     end
%             end
% 
%             if isempty(lowIndmid) || isempty(lowIndmid)
%                 error('error in crushing')
%             end
%             
%             if length(cracksXmid) >100
%                 pause(0.1)
%             end
%             
%             iceField.X = horzcat(iceField.X(1:lowIndmid), cracksXmid, iceField.X(highIndmid:end));
%             iceField.Y = horzcat(iceField.Y(1:lowIndmid), cracksYmid, iceField.Y(highIndmid:end));
%             iceField.x = horzcat(iceField.x(1:lowIndmid), cracksxmid, iceField.x(highIndmid:end));
%             iceField.y = horzcat(iceField.y(1:lowIndmid), cracksymid, iceField.y(highIndmid:end));           
%         end
%     end
% end