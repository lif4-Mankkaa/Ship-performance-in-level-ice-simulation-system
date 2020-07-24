function [ iceForce, iceField, nLoad ] = resolveContact( r, v, ship, iceField, reallyin, l_c, t, lin, dyna, cra, randomNumber, DB)
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
        cy3 = cos(r(3));
        sy3 = sin(r(3));
        
        [contacts, noOfContacts] = grouping(reallyin);

        %/////////////////////////////////////////////////////////////////////////
        %/////////////////////////////////////////////////////////////////////////
        %determination of intersection points for each contact
        
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
            %in case only one intersection point
            if xinte == xints & yinte == yints
                xinte = NaN;
                noOfContacts = noOfContacts - 1;
            end            
            
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
        %/////////////////////////////////////////////////////////////////////////
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
        
        %/////////////////////////////////////////////////////////////////////////
        %/////////////////////////////////////////////////////////////////////////
        %interpolation of the properties of the ice field for the contacts
        midX = r(1) + midx.*cy3 - midy.*sy3;
        midY = r(2) + midx.*sy3 + midy.*cy3;        
        %ice thickness at midpoints is interpolated from the iceField
%         hi = interp2(iceField.N, iceField.M, iceField.h, midX, midY, 'nearest');
%         sigmac = interp2(iceField.N, iceField.M, iceField.sigmac, midX, midY, 'nearest');
%         sigmaf = interp2(iceField.N, iceField.M, iceField.sigmaf, midX, midY, 'nearest');
        hi = interp1(iceField.N(100,:)-1, iceField.h(100,:), min(round(max(midX,0)),length(iceField.N(100,:))-1), 'nearest');
        sigmaf = interp1(iceField.N(100,:)-1, iceField.sigmaf(100,:), min(round(max(midX,0)),length(iceField.N(100,:))-1), 'nearest');
        sigmac = interp1(iceField.N(100,:)-1,iceField.sigmac(100,:), min(round(max(midX,0)),length(iceField.N(100,:))-1),'nearest');
        l_c = (iceField.E*hi.^3/(12*(1-0.3.^2)*1025*9.81)).^0.25;
% tic
%         hi = iceField.h(100,(min(round(max(midX,0)),length(iceField.N(100,:))-1))');
%         sigmaf = iceField.sigmaf(100,(min(round(max(midX,0)),length(iceField.N(100,:))-1))');
%         sigmac = iceField.sigmac(100,(min(round(max(midX,0)),length(iceField.N(100,:))-1))');
%         toc
        %/////////////////////////////////////////////////////////////////////////
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
        end
                
        %relative motion of the contact area midpoints. y(4:6) are the x-, y- and
        %rotational component of vcog. vrelt, vreln and creln1 correspond to the
        %like named forces in the paper.
        vrel = [v(1)-v(3).*midy(:) v(2)+v(3).*midx(:)];
        vrelt = dot(vrel, tangent,2);
        vreln = dot(vrel, normal,2);
        vreln1 = vreln .* cos(meanPhi);
        vreln2 = vreln .* sin(meanPhi);

        %/////////////////////////////////////////////////////////////////////////
        %/////////////////////////////////////////////////////////////////////////
        %determination of contact areas
        
        %initialization of area and meanPhi vectors
        area = NaN(noOfContacts,1);

        cmPhi = cos(meanPhi);
        smPhi = sin(meanPhi);
        tmPhi = tan(meanPhi);
        
        j = 1;
        
        LC_DB = -1*ones(noOfContacts,1);
        
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

            correctedDistance = (distance.*abs(tan(meanPhi(j)))>hi(j)).*(hi(j)/sin(meanPhi(j))) + (distance.*abs(tan(meanPhi(j)))<=hi(j)).*distance./abs(cos(meanPhi(j)));
            %distance to move ice nodes for area calculation
            moveDistance = correctedDistance-distance;
            %determination of corrected ice nodes for area calculation
            if sign(ship.y(contactswl_1(i))) ==sign(ship.y(contactswl_1(i+1)))            
                if sign(a(j))*(a(j)*xice0+b(j)*yice0+c(j)) == 0
                    xice = xice0 + moveDistance .* normal(j,1);
                    yice = yice0 + moveDistance .* normal(j,2);
                else            
                    xice = logical(sign(a(j)).*(a(j)*xice0+b(j)*yice0+c(j)) > 0) .* (xice0 - moveDistance .* normal(j,1)) + logical(sign(a(j)).*(a(j)*xice0+b(j)*yice0+c(j)) < 0) .* (xice0 + moveDistance .* normal(j,1));
                    yice = logical(sign(a(j)).*(a(j)*xice0+b(j)*yice0+c(j)) > 0) .* (yice0 - moveDistance .* normal(j,2)) + logical(sign(a(j)).*(a(j)*xice0+b(j)*yice0+c(j)) < 0) .* (yice0 + moveDistance .* normal(j,2));
                end
            else
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
            xin2 = vertcat(xint(i+1), xwl, xint(i));
            yin2 = vertcat(yint(i+1), ywl, yint(i));
            [ctrx, ctry] = areaCentroid(xint(i), yint(i), atan(a(j)*a0(j)), xin1, yin1, xin2, yin2);
            if isnan(ctrx)
                error('ctrx is empty')
            end
            %-------------------------------------------------------
            %crushing center and slope, crushing depth
            center_dist = (abs(a(j)*ctrx+b(j)*ctry+c(j))/sqrt(a(j)^2+b(j)^2));%distance from center to crushing line
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
            k1 = (yint(i)-yice0(1))/(xint(i)-xice0(1));
            k2 = (yint(i+1)-yice0(end))/(xint(i+1)-xice0(end));
            if abs(k1) == Inf && k2 ~= Inf
                xtip = xint(i);
                ytip = k2*(xtip - xint(i+1))+yint(i+1);
            elseif k1 ~= Inf && k2 == Inf
                xtip = xint(i+1);
                ytip = k1*(xtip - xint(i))+yint(i);
            elseif k1 == Inf && k2 == Inf
                xtip = nan;ytip = nan;
            else
                xtip = (k1*xint(i)-yint(i)+yint(i+1)-k2*xint(i+1))/(k1-k2);
                ytip = k1*(xtip - xint(i))+yint(i);
            end
                        
            aaa = ((xtip-xint(i))^2+(ytip-yint(i))^2)^0.5;
            bbb = ((xtip-xint(i+1))^2+(ytip-yint(i+1))^2)^0.5;
            ccc = ((xint(i)-xint(i+1))^2+(yint(i)-yint(i+1))^2)^0.5;
            the(j) = acosd((aaa^2+bbb^2-ccc^2)/(2*aaa*bbb));
            alp1 = acosd((aaa^2+ccc^2-bbb^2)/(2*aaa*ccc));
            alp2 = 180-the(j)-alp1;
            alp = max(alp1,alp2);
            posneg = (alp1<=alp2) + (alp1>alp2)*(-1);
            alpha_unknown = alp-(180-the(j))/2;
            alpha_unknown(isnan(alpha_unknown))=0;
            
%             alp1 = min(abs(alpha-alpha1),180-abs(alpha-alpha1));
%             alp2 = min(abs(alpha-alpha2),180-abs(alpha-alpha2));
%             alp = max(alp1,alp2);
%             the = abs(alpha1-alpha2);
%             if alp < the
%                 alp = 180 - alp;
%             end
%             alpha_unknown = 90-alp+0.5*the;
            k_DB(j) = tand(alpha_unknown)*posneg;
            
            %contact area
            area(j) = polyarea(xbound, ybound);
            
            %crushing depth
            if  (contactswl(i) ~= -10000) 
                if sign(ship.y(contactswl(i))) ~=sign(ship.y(contactswl(i+1)))
%                     depth1 = (abs(a(j)*66.1016+c(j))/sqrt(a(j)^2+b(j)^2));
%                     depth2 = ((xint(i)-xint(i+1))^2+(yint(i)-yint(i+1))^2)^0.5;
%                     depth(j) = (0.5*depth1*depth2/pi)^0.5;
                    depth(j) = (area(j)*2/pi)^0.5;
                    LC_DB(j) = 0.33;
                else
                    depth(j) = (abs(a(j)*ctrx+b(j)*ctry+c(j))/sqrt(a(j)^2+b(j)^2))/LC_DB(j);
                end
            else
                depth(j) = (abs(a(j)*ctrx+b(j)*ctry+c(j))/sqrt(a(j)^2+b(j)^2))/LC_DB(j);
            end
%             
%             scatter(xice0,yice0);hold on
%             scatter(xint(i:i+1),yint(i:i+1))
%             plot(iceField.x,iceField.y,ship.x,ship.y)
%             xlim([min(xint(i:i+1))-0.1,max(xint(i:i+1))+0.1])
%             ylim(([min(yint(i:i+1))-0.1,max(yint(i:i+1))+0.1]))
%             hold off            
            %-------------------------------------------------------
            
            midX(j) = r(1) + ctrx.*cy3 - ctry.*sy3;
            midY(j) = r(2) + ctrx.*sy3 + ctry.*cy3;
                        
            if isnan(midX(j))
                disp('ri le gou')%error
            end            
            
            if sign(ship.y(contactswl_1(i))) ~= sign(ship.y(contactswl_1(i+1)))
                area(j) = area(j)/meanPhi(j);
            end            
            
            j = j + 1;
            
        end
        
        %/////////////////////////////////////////////////////////////////////////
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

        nice = length(iceField.X);%number of ice nodes
        distL = 0.05;%discretization length of ice, HOIDA SITEN ETT?SÄÄTYY KAIKKEEN KÄYTTÖÖNYHDEST?PAIKASTA
        dang = distL./R;%incrimental angle (absolute value)
        cracks = NaN*ones(ceil(2*pi/min(dang)), 2, length(R));%initialize matrix for X and Y coordinates of crack points
        crStEnd = zeros(length(R), 2);
        bendCounter = 1;
        
        wedgeAng = nan(length(R),1);
        
        B_ref = [0.005:0.001:0.04,0.042:0.002:0.05,0.055:0.005:0.1];
        E = iceField.E;
        SF = sigmac/1e6.*(E/1e9./hi).^0.5;
        SF_H = sigmac/1e6.*(E/1e9./hi).^0.25;
        
        for i = 1:length(R)
            
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
            angleInd = round((meanPhi(i)/pi*180-20)/5)+1;
%             %Tricky
%             if angleInd == 14
%                 angleInd = 13;
%             end
                        %
            DBcolumeInd = 1;
            if theta>pi/3 && theta<5/6*pi
                  thetaInd = round((theta-pi/3)/(pi/60)+1);
                  if angleInd >13
                      rbreak = Inf;
                      Edge1 = 1;
                      Edge2 = 1;
                      center = 1;
                  else
                      rbreak_ind = zeros(100,1);
                      jj = 0;
                      while length(rbreak_ind)==100 || DBcolumeInd==0
                            db = dbdb.DB0;
%                           switch round(abs(10*k_DB(i)))-jj
%                             case 0
%                                 db = dbdb.DB0;
%                             case 1
%                                 db = dbdb.DB1;
%                             case 2
%                                 db = dbdb.DB2;
%                             case 3
%                                 db = dbdb.DB3;
%                             case 4
%                                 db = dbdb.DB4;
%                             case 5
%                                 db = dbdb.DB5;
%                             case 6
%                                 db = dbdb.DB6;
%                             case 7
%                                 db = dbdb.DB7;
%                             case 8
%                                 db = dbdb.DB8;
%                             case 9
%                                 db = dbdb.DB9;
%                             otherwise
%                                 db = dbdb.DB10;
%                           end
                          
%                           %tricky
%                           if angleInd == 13
%                               db = dbdb.DB0;
%                           end
                          
                        jj = jj+1;
                        if jj >2
                            jj
                        end
                        
%                         kk = 0;
%                         DBcolumeInd = nan;
%                         while isnan(DBcolumeInd)
%                             thetaInd = thetaInd-kk;
                        comparedata = db(thetaInd).stress(angleInd,:);
                        comparedataH = db(thetaInd).stressH(angleInd,:);
                        rbreak_ind = find(comparedata*SF(i)+comparedataH*SF_H(i)<sigmaf(i));
                        comparedata(isnan(comparedata))=[];
                        DBcolumeInd = min(length(rbreak_ind)+1,length(comparedata));
%                             kk = kk+1;
%                             disp('thetaIndtrick')
%                         end
                      end
                        if length(rbreak_ind) == length(comparedata)
                            rbreak = Inf;
                            Edge1 = 1;Edge2 = 1;center = 1;
                        else
                            rbreak = B_ref(length(rbreak_ind)+1)*l_c(i);  
                        
                                                                    
                            if sign(k_DB(i)) == 1
                                Edge1 = db(thetaInd).edge1(angleInd,DBcolumeInd)*l_c(i);
                                Edge2 = db(thetaInd).edge2(angleInd,DBcolumeInd)*l_c(i);
                            else
                                Edge1 = db(thetaInd).edge2(angleInd,DBcolumeInd)*l_c(i);
                                Edge2 = db(thetaInd).edge1(angleInd,DBcolumeInd)*l_c(i);
                            end
                            center = db(thetaInd).center(angleInd,DBcolumeInd)*l_c(i);
                        end
%                         %randomize cusp geometry
%                         Edge1 = 0.27/0.94*randomNumber*Edge1+Edge1;
%                         Edge2 = 0.27/0.94*randomNumber*Edge2+Edge2;
%                         center = 0.27/0.94*randomNumber*center+center;
                  end   
                    
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
                    
                    theta = N*dang(i);
                    coef = [0 0 1;0.25*theta^2 0.5*theta 1;theta^2 theta 1]\[Edge1;center;Edge2];
                    
                    cm = tril(ones(N));%coefficient matrix
                    th = lowAng+cm*(dang(i)*ones(N,1));
                    rrr = coef(1)*(th-lowAng).^2+coef(2)*(th-lowAng)+coef(3);

                    cracksx = midX(i)+rrr.*cos(th);
                    cracksy = midY(i)+rrr.*sin(th);
             else
                  db = dbdb.DB0;
                  thetaInd = 1;
                  if angleInd >13
                      rbreak = Inf;
                      rrr = 1;
                  else 
                    rbreak_ind = find(db(thetaInd).stress(angleInd,:)*SF(i)<sigmaf(i));
                    
                    if length(rbreak_ind) == 51
                        rbreak = Inf;
                    else
                        rbreak = B_ref(length(rbreak_ind)+1)*l_c(i); 
                    end
                    DBcolumeInd = min(length(rbreak_ind)+1,51);

                    rrr = db(thetaInd).center(angleInd,DBcolumeInd)*l_c(i);
%                     %randomize cusp geometry
%                     rrr = 0.27/0.94*randomNumber*rrr+rrr;
                  end 
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
            if depth(i) >= rbreak
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
            else
                cracks(:,:,bendCounter) = NaN;
                R(i) = NaN;
                Redge(i) = NaN;
                crStEnd(i) = NaN;
            end
          end

        
        retain = ~isnan(R);
        RR = R;
        if ~isempty(R)
                    midX = midX(retain);
                    midY = midY(retain);
        end
        R = R(retain);
        crStEnd = crStEnd(retain,:);
        cracks = cracks(:,:,retain);
        
         %crushing at midship area
        iceField.x = 1/(cos(r(3))^2+sin(r(3))^2)*(iceField.X*cos(r(3))+iceField.Y*sin(r(3)))-(r(1)*cos(r(3))+r(2)*sin(r(3)));
        iceField.y = 1/(cos(r(3))^2+sin(r(3))^2)*(-iceField.X*sin(r(3))+iceField.Y*cos(r(3)))+(r(1)*sin(r(3))-r(2)*cos(r(3)));
        
         for i = 1:length(xint)/2
             if abs(xint(2*i-1)-xint(2*i))>0.00001
                xintsort = sort([xint(2*i-1),xint(2*i)]);
                if (xintsort(1)>=ship.mids && xintsort(2)<=ship.mide)% || (xintsort(2)>=ship.mids && xintsort(2)<=ship.mide)
                    n = max(ceil(abs(xint(2*i)-xint(2*i-1))/0.1),1);%number of nodes of new crack
                    if diff(xintsort)>0.1
                        cracksxmid = xint(2*i-1):(xint(2*i)-xint(2*i-1))/n:xint(2*i);
                        cracksymid = yint(2*i-1)*ones(1,n+1);
                    else
                        cracksxmid = xint(2*i);
                        cracksymid = yint(2*i);
                    end
                    cracksXmid = r(1) + cracksxmid*cy3 - cracksymid*sy3;
                    cracksYmid = r(2) + cracksxmid*sy3 + cracksymid*cy3;
                    if yint(2*i-1)<0
                        lowIndmidall = find(iceField.x(1:end-4)<=xint(2*i-1)&iceField.y(1:end-4)<yint(2*i-1));
                        if isempty(lowIndmidall)
                            closest = min(iceField.x(iceField.y<yint(2*i-1)) - xint(2*i-1))+xint(2*i-1);
                            lowIndmidall = find(iceField.x == closest);
                        end
                        lowIndmid = lowIndmidall(end);
                        highIndmid = find(iceField.x(1:end-4)>xint(2*i)&iceField.y(1:end-4)<yint(2*i-1),1);
                    else
                        lowIndmidall = find(iceField.x(1:end-4)>xint(2*i-1)&iceField.y(1:end-4)>yint(2*i-1));                      
                        lowIndmid = lowIndmidall(end);
                        highIndmid = find(iceField.x(1:end-4)<xint(2*i)&iceField.y(1:end-4)>yint(2*i-1),1);
                        if isempty(highIndmid)
                            closest = min(iceField.x(iceField.y>yint(2*i-1)) - xint(2*i-1))+xint(2*i-1);
                            highIndmid1 = find(iceField.x == closest);
                            highIndmid = highIndmid1(1);
                        end                    
                    end
                        iceField.X = horzcat(iceField.X(1:lowIndmid), cracksXmid, iceField.X(highIndmid:end));
                        iceField.Y = horzcat(iceField.Y(1:lowIndmid), cracksYmid, iceField.Y(highIndmid:end));
                        iceField.x = horzcat(iceField.x(1:lowIndmid), cracksxmid, iceField.x(highIndmid:end));
                        iceField.y = horzcat(iceField.y(1:lowIndmid), cracksymid, iceField.y(highIndmid:end));
                end    
             end
         end
%         end
        
        %/////////////////////////////////////////////////////////////////////////
        %/////////////////////////////////////////////////////////////////////////
        %forces in the x- and y direction
%         if jjj ~= 1
        xi = repmat([1 0], length(fH), 1);
        yi = repmat([0 1], length(fH), 1);
        Fx = dot(tangent,xi,2).*fH + dot(normal,xi,2).*FH;
        Fy = dot(tangent,yi,2).*fH + dot(normal,yi,2).*FH;
        Mz = Fx .* midy(:) + Fy .* midx(:);
        %compute total x and y direction forces and moment
        iceForce = [sum(Fx(isfinite(Fx))); sum(Fy(isfinite(Fy))); sum(Mz(isfinite(Mz)))];    
        nLoad = [midx(:) midy(:) cLength(:) Fcr(:) FH(:) FV(:) fH(:) fV(:) LC_DB(:) sigmac(:) sigmaf(:) hi(:) RR(:) Redge(:) wedgeAng(:) k_DB(:)];
%         else
%             Fx = 0;
%             Fy = 0;
%             Mz = 0;
%             iceForce = [0;0;0]
%             nLoad = [0,0,0,0,0,0,0,0,0,0,0,0,0,0]
%         end
%          end
%      else
%          iceForce = [0;0;0];
%          nLoad = [];
%      end
end