
function [coe,r_vec_SOI_Entry_Ia]  = func_copy_MAIN_Planetary_Approach(v_inf_I,zp,i,Primary)
% coe = (e,a,i,argp,RAAN,TA,h)

% r_entry  = r_vec_SOI_Entry_Ia
% clear;close all
% clc;
% % INPUTS
% Primary.Name = 'Earth';   % Planet
% v_inf_I = [-3;5;2];      %[km/s] Excess Velocity Vector at SOI Entrance in Inertial, Planetocentric Coordinates
% zp = 10000;               %[km] Altitude at periapsis
% i = 30; % mod. 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ORBIT DETERMINATION
% Planet's Properties
planetparam = func_Primary_Parameters_RevC(Primary);
RE = planetparam.RE;
RP = planetparam.RP;
mu = planetparam.mu;
r_SOI = planetparam.r_SOI;
% [RE,RP,mu,~,~,~,~,~,r_SOI] = func_Primary_Parameters_RevB(Primary); % Don't Edit

% if strcmp(Primary.Name,'Venus')||strcmp(Primary.Name,'Mars')||strcmp(Primary.Name,'Jupiter')||strcmp(Primary.Name,'Saturn')
%     ColorTraj = [0,.7,1];
%     ColorPeri = [1,1,0];
% else
%     ColorTraj = [1,1,0];
%     ColorPeri = [0,.7,1];
% end

v1 = v_inf_I(1);
v2 = v_inf_I(2);
v3 = v_inf_I(3);

v = norm(v_inf_I);
i_min = asind(abs(v3)/v);
i_max = 180-asind(abs(v3)/v);

% i = i_min:(i_max-i_min)/10:i_max;
fprintf('Orbital Inclination %.1f <= i <= %.1f deg\n',i_min,i_max)

for k=1:size(i,2)
    if i(k)<i_min || i(k)>i_max 
        fprintf('Choose a valid orbital inclination.\n\n')
        return
    end
    rp = zp + RE;             %[km] Periapsis Radius
    a = mu/v^2;               %[km] Semi-major Axis
    e = rp/a+1;               %[--] Eccentricity
    h = mu/v*sqrt(e^2-1);     %[km^2/s] Angular Momentum Vector Norm
    En = v^2/2;               %[km^2/s^2] Orbital Energy
    Turn_Angle = 2*asind(1/e);%[deg] Turn Angle
    Aiming_Radius = a*sqrt(e^2-1); %[km] Aiming Radius
    TA_inf_Entry = -acosd(-1/e);    %[deg]
    TA_inf_Exit  =  acosd(-1/e);    %[deg]

    % Angular Momentum Vector
    h3(k)  = h*cosd(i(k));
    
    if v1~=0
        discr = v^2*sind(i(k))^2-v3^2;
        if abs(v^2*sind(i(k))^2-v3^2) < 1E-12
            discr = 0;
        end
        h2a(k) = h*(-v2*v3*cosd(i(k)) + v1*sqrt(discr))/(v1^2+v2^2);
        h1a(k) = -v2/v1*h2a(k) - v3/v1*h3(k);
    
        h2b(k) = h*(-v2*v3*cosd(i(k)) - v1*sqrt(discr))/(v1^2+v2^2);
        h1b(k) = -v2/v1*h2b(k) - v3/v1*h3(k);
     
        h_vec_a_I(:,k) = [h1a(k);h2a(k);h3(k)];
        h_vec_b_I(:,k) = [h1b(k);h2b(k);h3(k)];
    end
    if v2~=0
        discr = v^2*sind(i(k))^2-v3^2;
        if abs(v^2*sind(i(k))^2-v3^2) < 1E-12
            discr = 0;
        end
        h1a(k) = h*(-v1*v3*cosd(i(k)) - v2*sqrt(discr))/(v1^2+v2^2);
        h2a(k) = -v1/v2*h1a(k) - v3/v2*h3(k);
    
        h1b(k) = h*(-v1*v3*cosd(i(k)) + v2*sqrt(discr))/(v1^2+v2^2);
        h2b(k) = -v1/v2*h1b(k) - v3/v2*h3(k);
     
        h_vec_a_I(:,k) = [h1a(k);h2a(k);h3(k)];
        h_vec_b_I(:,k) = [h1b(k);h2b(k);h3(k)];
    end
    
    if v1==0 && v2==0
        i = 90;
        h3 = 0;
    end
    
    % p3-versor (Perifocal RF)
    p3a_I(:,k) = h_vec_a_I(:,k)/h;
    p3b_I(:,k) = h_vec_b_I(:,k)/h;
    
    % Line of Nodes
    N_vec_a_I(:,k) = [-h2a(k)/(h*sind(i(k)));h1a(k)/(h*sind(i(k)));0];
    N_vec_b_I(:,k) = [-h2b(k)/(h*sind(i(k)));h1b(k)/(h*sind(i(k)));0];
    
    % Right Ascension of the Ascending Node (RAAN)
    if N_vec_a_I(2,k)>0
        RAANa(k) = acosd(-h2a(k)/(h*sind(i(k))));
    else
        RAANa(k) = 360-acosd(-h2a(k)/(h*sind(i(k))));
    end
    if N_vec_b_I(2,k)>0
        RAANb(k) = acosd(-h2b(k)/(h*sind(i(k))));
    else
        RAANb(k) = 360-acosd(-h2b(k)/(h*sind(i(k))));
    end
    
    % Eccentricity Vector
    e_vec_a_I(:,k) = [(v2*h3(k) -v3*h2a(k))/mu+v1/v
                      (v3*h1a(k)-v1*h3(k) )/mu+v2/v
                      (v1*h2a(k)-v2*h1a(k))/mu+v3/v];
   
    e_vec_b_I(:,k) = [(v2*h3(k) -v3*h2b(k))/mu+v1/v
                      (v3*h1b(k)-v1*h3(k) )/mu+v2/v
                      (v1*h2b(k)-v2*h1b(k))/mu+v3/v];
    
    % p1-versor (Perifocal RF)
    p1a_I(:,k) = e_vec_a_I(:,k)/e;
    p1b_I(:,k) = e_vec_b_I(:,k)/e;
    
    % p2-versor (Perifocal RF)
    p2a_I(:,k) = cross(p3a_I(:,k),p1a_I(:,k));
    p2b_I(:,k) = cross(p3b_I(:,k),p1b_I(:,k));
    
    % Argument of Periapsis
    if e_vec_a_I(3,k)>0
        ARGPa(k) = acosd(N_vec_a_I(:,k)'*p1a_I(:,k));
    else
        ARGPa(k) = 360-acosd(N_vec_a_I(:,k)'*p1a_I(:,k));
    end
    if e_vec_b_I(3,k)>0
        ARGPb(k) = acosd(N_vec_b_I(:,k)'*p1b_I(:,k));
    else
        ARGPb(k) = 360-acosd(N_vec_b_I(:,k)'*p1b_I(:,k));
    end
    
    % Coordinate Transformation Matrix between PCI and P-RF
    R_PaI(:,:,k) = R3(ARGPa(k))*R1(i(k))*R3(RAANa(k));
    R_PbI(:,:,k) = R3(ARGPb(k))*R1(i(k))*R3(RAANb(k));
    
    % Eccentricity on P-RF Coordinates
    e_vec_a_P(:,k) = R_PaI(:,:,k)*e_vec_a_I(:,k);
    e_vec_b_P(:,k) = R_PaI(:,:,k)*e_vec_b_I(:,k);

    % V_inf in P-RF Coordinates
    v_inf_Pa(:,k)  = R_PaI(:,:,k)*v_inf_I;
    v_inf_Pb(:,k)  = R_PbI(:,:,k)*v_inf_I;
    
    % Ang. Mom. Vector in P-RF Coordinates
    h_vec_Pa(:,k)  = R_PaI(:,:,k)*h_vec_a_I(:,k);
    h_vec_Pb(:,k)  = R_PbI(:,:,k)*h_vec_b_I(:,k);
    
    % SOI Entrance
    TA_SOI_Entry = -acosd(1/e*(h^2/mu/r_SOI -1));
    r_vec_SOI_Entry_P = r_SOI*[cosd(TA_SOI_Entry); sind(TA_SOI_Entry); 0];   
    r_vec_SOI_Entry_Ia(:,k) = R_PaI(:,:,k)'*r_vec_SOI_Entry_P;
    r_vec_SOI_Entry_Ib(:,k) = R_PbI(:,:,k)'*r_vec_SOI_Entry_P;

    DEC_Entry_a(k) =  90 - acosd(dot(r_vec_SOI_Entry_Ia(:,k)/r_SOI,[0;0;1]));
    RA_Entry_a(k)   = atan2d(r_vec_SOI_Entry_Ia(2,k),r_vec_SOI_Entry_Ia(1,k));
    if RA_Entry_a(k)<0
        RA_Entry_a(k)=RA_Entry_a(k)+360;
    end
    DEC_Entry_b(k) =  90 - acosd(dot(r_vec_SOI_Entry_Ib(:,k)/r_SOI,[0;0;1]));
    RA_Entry_b(k)   = atan2d(r_vec_SOI_Entry_Ib(2,k),r_vec_SOI_Entry_Ib(1,k));
    if RA_Entry_b(k)<0
        RA_Entry_b(k)=RA_Entry_b(k)+360;
    end

    % SOI Exit
    TA_SOI_Exit  =  acosd(1/e*(h^2/mu/r_SOI -1));
    r_vec_SOI_Exit_P = r_SOI*[cosd(TA_SOI_Exit); sind(TA_SOI_Exit); 0];
    r_vec_SOI_Exit_Ia(:,k) = R_PaI(:,:,k)'*r_vec_SOI_Exit_P;
    r_vec_SOI_Exit_Ib(:,k) = R_PbI(:,:,k)'*r_vec_SOI_Exit_P;

    DEC_Exit_a(k) =  90 - acosd(dot(r_vec_SOI_Exit_Ia(:,k)/r_SOI,[0;0;1]));
    RA_Exit_a(k)   = atan2d(r_vec_SOI_Exit_Ia(2,k),r_vec_SOI_Exit_Ia(1,k));
    if RA_Exit_a(k)<0
        RA_Exit_a(k)=RA_Exit_a(k)+360;
    end
    DEC_Exit_b(k) =  90 - acosd(dot(r_vec_SOI_Exit_Ib(:,k)/r_SOI,[0;0;1]));
    RA_Exit_b(k)   = atan2d(r_vec_SOI_Exit_Ib(2,k),r_vec_SOI_Exit_Ib(1,k));
    if RA_Exit_b(k)<0
        RA_Exit_b(k)=RA_Exit_b(k)+360;
    end
% func_OEtoRV(e,a,i,argp,RAAN,TA,mu)
coe = [e, a, i(k), ARGPa, RAANa,TA_SOI_Entry, h];
% coe2 = [e, a, i(k), ARGPb, RAANb,TA_SOI_Entry];
    fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
    fprintf('   Orbital Inclination:     i = %.1f deg\n',i(k))
    fprintf(' Altitude of Periapsis:    zp = %.1f km\n',zp)
    fprintf('          Excess Speed: v_inf = %.3f km/s\n',v)
    fprintf(' Characteristic Energy:    C3 = %.1f km^2/s^2\n',v^2)
    fprintf('        Orbital Energy:     E = %.1f km^2/s^2\n',En)
    fprintf('       Semi-Major Axis:     a = %.1f km\n',a)
    fprintf('          Eccentricity:     e = %.3f\n',e)
    fprintf(' Angular Momentum Norm:     h = %.1f km^2/s\n',h)
    fprintf('            Turn Angle:     d = %.1f deg\n',Turn_Angle)
    fprintf('         Aiming Radius:     D = %.1f km\n',Aiming_Radius)
    fprintf('Right Asc.of Asc. Node:  RAAN = %.1f deg  or  %.1f deg\n',RAANa(k),RAANb(k))
    fprintf(' Argument of Periapsis:  ARGP = %.1f deg  or  %.1f deg\n\n',ARGPa(k),ARGPb(k))
    fprintf('        @ SOI Entrance:    RA = %.2f deg  or  %.2f deg\n',RA_Entry_a(k),RA_Entry_b(k))
    fprintf('                          DEC = %.2f deg  or  %.2f deg\n',DEC_Entry_a(k),DEC_Entry_b(k))
    fprintf('        @ SOI     Exit:    RA = %.2f deg  or  %.2f deg\n',RA_Exit_a(k),RA_Exit_b(k))
    fprintf('                          DEC = %.2f deg  or  %.2f deg\n',DEC_Exit_a(k),DEC_Exit_b(k))
    fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')

    TA = TA_SOI_Entry:(TA_SOI_Exit-TA_SOI_Entry)/1000:TA_SOI_Exit;

    for j=1:size(TA,2)
        r(j) = h^2/mu/(1+e*cosd(TA(j)));
        r_P(:,j) = r(j)*[cosd(TA(j)); sind(TA(j)); 0];
        r_Ia(:,j,k) = R_PaI(:,:,k)'*r_P(:,j);
        r_Ib(:,j,k) = R_PbI(:,:,k)'*r_P(:,j);
    end
    r_PER_P  = [h^2/mu/(1+e);0;0];
    r_PER_Ia(:,k) = R_PaI(:,:,k)'*r_PER_P;
    r_PER_Ib(:,k) = R_PbI(:,:,k)'*r_PER_P;
end
end


% Color_Entry = [0,0.7,0];
% Color_Exit  = [1,0.4,0];
% 
% figure(1),hold on
% set(gcf,'units','normalized','position',[0.0,0.0,1.0,0.9])
% 
% subplot(1,2,1),hold on
% tlt1 = title('Sphere of Influence (SOI) & Celestial Equatorial Plane');
% set(tlt1,'Color','w')
% set(gca,'Color','k')
% set(gcf,'Color','k')    
% set(gca,'XTick',[], 'YTick', [], 'ZTick', [])
% ax = gca;
% ax.XColor = 'none'; ax.YColor = 'none'; ax.ZColor = 'none';
% 
% axis('equal')
% view([-7,16])
% [xSOI,ySOI,zSOI] = ellipsoid(0,0,0,r_SOI,r_SOI,r_SOI,60);
% SOI_Plot = surface(xSOI,ySOI,zSOI);
% SOI_Plot.FaceColor = 'b';
% SOI_Plot.FaceAlpha = 0.15;
% SOI_Plot.EdgeColor = 'none';
% xlim([-r_SOI,r_SOI])
% ylim([-r_SOI,r_SOI])
% zlim([-r_SOI,r_SOI])
% 
% phi = 0:6:360;
% for j=1:size(phi,2)
%     A_I(:,j) = r_SOI*[cosd(phi(j)),sind(phi(j)),0];
% end
% Eq_Plane = patch(A_I(1,:),A_I(2,:),A_I(3,:),[1,0,1]); 
% Eq_Plane.FaceColor = [0,0.5,1];
% Eq_Plane.EdgeColor = 'none';
% Eq_Plane.FaceAlpha = 0.2; 
% 
% 
% for k=1:size(i,2)
%     plot3(r_Ia(1,:,k),r_Ia(2,:,k),r_Ia(3,:,k),'color',ColorTraj)
%     plot3(r_Ib(1,:,k),r_Ib(2,:,k),r_Ib(3,:,k),'color',ColorTraj)
% end
%
% plot3(r_vec_SOI_Entry_Ia(1,:),r_vec_SOI_Entry_Ia(2,:),r_vec_SOI_Entry_Ia(3,:),'marker','s','markeredgecolor',Color_Entry,'markerfacecolor',Color_Entry,'markersize',6,'linestyle','none');
% plot3(r_vec_SOI_Entry_Ib(1,:),r_vec_SOI_Entry_Ib(2,:),r_vec_SOI_Entry_Ib(3,:),'marker','s','markeredgecolor',Color_Entry,'markerfacecolor',Color_Entry,'markersize',6,'linestyle','none');
% 
% plot3(r_vec_SOI_Exit_Ia(1,:),r_vec_SOI_Exit_Ia(2,:),r_vec_SOI_Exit_Ia(3,:),'marker','s','markeredgecolor',Color_Exit,'markerfacecolor',Color_Exit,'markersize',6,'linestyle','none');
% plot3(r_vec_SOI_Exit_Ib(1,:),r_vec_SOI_Exit_Ib(2,:),r_vec_SOI_Exit_Ib(3,:),'marker','s','markeredgecolor',Color_Exit,'markerfacecolor',Color_Exit,'markersize',6,'linestyle','none');
% 
% subplot(1,2,2),hold on
% tlt2 = title('Zoom-In on Planet with Locus of Periapsis');
% set(tlt2,'Color','w')
% set(gca,'Color','k')
% set(gcf,'Color','k')
% ax = gca;
% ax.XColor = 'w'; ax.YColor = 'w'; ax.ZColor = 'w';
% 
% axis('equal')
% [xPRI,yPRI,zPRI] = ellipsoid(0,0,0,RE/RE,RE/RE,RP/RE,60);
% PRI_Plot = surface(xPRI,yPRI,-zPRI);
% PRI_Plot.FaceColor = 'texturemap';
% PRI_Plot.FaceAlpha = 1.0;
% PRI_Plot.EdgeColor = 'none';
% PRI_Plot.CData     = Primary.Figure;
% xlim([-2*rp/RE,2*rp/RE])
% ylim([-2*rp/RE,2*rp/RE])
% zlim([-2*rp/RE,2*rp/RE])
% view([0,0])
% 
% xlabel('X - Planet''s Radii')
% ylabel('Y - Planet''s Radii')
% zlabel('Z - Planet''s Radii')
% for k=1:size(i,2)
%     plot3(r_Ia(1,:,k)/RE,r_Ia(2,:,k)/RE,r_Ia(3,:,k)/RE,'color',ColorTraj)
%     plot3(r_Ib(1,:,k)/RE,r_Ib(2,:,k)/RE,r_Ib(3,:,k)/RE,'color',ColorTraj)
% end
% plot3(r_PER_Ia(1,:)/RE,r_PER_Ia(2,:)/RE,r_PER_Ia(3,:)/RE,'marker','s','markeredgecolor',ColorPeri,'markerfacecolor',ColorPeri,'markersize',5,'linestyle','none');
% plot3(r_PER_Ib(1,:)/RE,r_PER_Ib(2,:)/RE,r_PER_Ib(3,:)/RE,'marker','s','markeredgecolor',ColorPeri,'markerfacecolor',ColorPeri,'markersize',5,'linestyle','none');
% 
% 
% %% Extra Functions
% function [RE,RP,mu,TE,J2,Color_Primary,Color_Primary_Neg,Map_Primary,r_SOI] = func_Primary_Parameters1(Primary)
%     % Planetary Parameters
%     % RE: Equatorial Radius                [km]
%     % RP: Polar Radius                     [km]
%     % f:  flatness                         [-]
%     % mu: Gravitational Parameter          [km^3/s^2]
%     % TE: Sidereal Period                  [s]
%     % J2: Second Zonal Spherical Harmonic  [-]
%     if strcmp(Primary,'Mercury')
%         RE = 2440;
%         mu = 22032;
%         TE = 58.65*86400;
%         J2 = 60E-6;
%         f  = 0;
%         Color_Primary = [.55,.55,.56];
%         r = 57.91e6;
%     elseif strcmp(Primary,'Venus')
%         RE = 6052;
%         mu = 324859;
%         TE = -243*86400;
%         J2 = 4.458E-6;
%         f = 0;
%         Color_Primary = [1,.5,0];
%         r = 108.2e6;
%     elseif strcmp(Primary,'Earth')    
%         RE = 6378.1363;
%         mu = 398600.4415;
%         TE = 0.99726968*86400;
%         J2 = 1.0826269E-3;
%         f = 0.0033528131;
%         Color_Primary = [.2,.2,1];
%         r = 149.6E6;
%     elseif strcmp(Primary,'Mars')   
%         RE = 3396;
%         mu = 42828;
%         TE = 24.62*3600;
%         J2 = 1.96045E-3;
%         f = 0.00648;
%         Color_Primary = [1,.4,.1];
%         r = 227.9e6;
%     elseif strcmp(Primary,'Jupiter')
%         RE = 71490;
%         mu = 126686534;
%         TE = 9.925*3600;
%         J2 = 14.736E-3;
%         f = 0.06487;
%         Color_Primary = [.9,.7,.2];
%         r = 778.6e6;
%     elseif strcmp(Primary,'Saturn') 
%         RE = 60270;
%         mu = 37931187;
%         TE = 10.66*3600;
%         J2 = 16.298E-3;
%         f = 0.09796;
%         Color_Primary = [1,.9,.6];
%         r = 1433e6;
%     elseif strcmp(Primary,'Uranus') 
%         RE = 25560;
%         mu = 5793939;
%         TE = -17.24*3600;
%         J2 = 3.34343E-3;
%         f = 0.02293;
%         Color_Primary = [0,.9,1]; 
%         r = 2872e6;
%     elseif strcmp(Primary,'Neptune')
%         RE = 24764;
%         mu = 6836529;
%         TE = 16.21*3600;
%         J2 = 3.411E-3;
%         f = 0.01708;
%         Color_Primary = [0,.5,1];
%         r = 4495e6; 
%     elseif strcmp(Primary,'Moon')
%         RE = 1737;
%         mu = 4905;
%         TE = 27.32*86400;
%         J2 = 202.7E-6;
%         f = 0.0012;
%         Color_Primary = [.6,.61,.6];
%         r = 384.4e3; 
%     end
%         Primary
%         % Polar Radius
%         RP = RE*(1-f);
% 
%         % Colors
%         HSV_CP = rgb2hsv(Color_Primary);
%         if HSV_CP(2)<0.5
%             if HSV_CP(1)<0.5
%                 Color_Primary_Neg = hsv2rgb([HSV_CP(1)+0.5,1,1]);
%             else
%                 Color_Primary_Neg = hsv2rgb([HSV_CP(1)-0.5,1,1]);
%             end
%         else
%             Color_Primary_Neg = [1,1,1]-Color_Primary;
%         end
% 
%         % Planet's map
%         image_file  = sprintf('Map_%s.jpg',Primary);% Changed the filename Convenction
%         Map_Primary = imread(image_file);           % now "Map_planet"
% 
%         % Radius of Sphere of Influence (SOI)
%         mu_Sun = 1.3271244e11; %[km^3/s^2]
%         r_SOI = r*(mu/mu_Sun)^(2/5); %[km]
% end


function R1 = R1(a)
    R1 = [1,      0,       0;
          0, cosd(a),sind(a);
          0,-sind(a),cosd(a)];
end

function R2 = R2(a)
    R2 = [cosd(a),0,-sind(a);
                0,1,       0;
          sind(a),0, cosd(a)];
end

function R3 = R3(a)
    R3 = [ cosd(a),sind(a),0;
          -sind(a),cosd(a),0;
                 0,      0,1];
end