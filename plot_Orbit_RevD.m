
function [oiplot, oilplot,oiplot2, oilplot2, r_new] = plot_Orbit_RevD(r,v,Planet,UIaxishand,UIaxishand2)
% Plots the orbit based on r and v. 
% Plots to a figure, then outputs the figure handle.
% Theta inf is inputting in deg. 
% Code is formerly rads.
% r_new is in ECI coordinates

fxname = "plot_Orbit_RevC.m";
fprintf('\n\nStart "%s" \n',fxname)
planetparam= func_Primary_Parameters_RevC(Planet);
mu = planetparam.mu;

coe                = coe_from_sv(r,v,mu,false); % "false" suppresses output
h = coe(1); 
e = coe(2);
orbital_energy     = -mu^2/2/h^2*(1-e^2);
fprintf('Specific Orbital Energy,eps = %3.4f km^2/s^2\n',orbital_energy) 

%   coe  - vector of orbital elements [h e RA incl w TA a] - angle in rad.
% Euler Angles; RAAN (Omega), inclination (i), arguement of perigee (w)
% With the Euler Angles, make the rotation matrix.
% Sequence alpha beta gamma = omega, inclination, argument of perigree.

% All points in orbit are in flat plane (r,theta (true anonmaly))

if e>1 % Hyperbolic Orbit
    %theta inf = 
    thetainf_trueAnonofAsymtote = 180-acosd(1/e); %Vallado 2-29 5E
    fprintf('True anomaly of Asymptote, Theta_inf: %.0f deg \n',thetainf_trueAnonofAsymtote) %*180/pi
     k= 0.995; % scaling
    tru_anon_range = linspace(-thetainf_trueAnonofAsymtote*k,thetainf_trueAnonofAsymtote*k,100);
    tru_anon_range(end+1) = coe(6)*180/pi; % Adding True Anomaly to the end.
    rp =coe(7)*(e-1);
    r_orb = (h^2/mu) ./ ( 1 + e*cosd(tru_anon_range));
   
    x_orb = r_orb.*cosd(tru_anon_range);
    y_orb = r_orb.*sind(tru_anon_range);
else    
    tru_anon_range = 0:pi/180:2*pi; % rad.
    tru_anon_range(end+1) = coe(6); % Adding True Anomaly to the end.
    r_orb = h^2./mu *(1 ./ ( 1 + e*cos(tru_anon_range)));
    % Conver to rectangular coordinates. 
    x_orb = r_orb.*cos(tru_anon_range);
    y_orb = r_orb.*sin(tru_anon_range);
    fprintf('Plotting Non-Hyperbolic Orbit \n')
end

z_orb = zeros(size(x_orb));

% Form the Rotation Matrix.
RA   = coe(3); % Omega_RAAN, rads 
incl = coe(4); % Inclination, rad
w    = coe(5); % omega_Arguement of Perigee.

% New Rotation Matrix
%...Equation 4.34:% R3 (phi)
R3_W = [ cos(RA)  sin(RA)  0
        -sin(RA)  cos(RA)  0
            0       0      1 ];
R1_i = [    1       0      0
            0   cos(incl)  sin(incl)
            0  -sin(incl)  cos(incl)];
R3_w = [ cos(w)  sin(w)    0 
        -sin(w)  cos(w)    0
           0       0       1 ];
rot_to_PF  = R3_w*R1_i*R3_W; % Coord. Transf. Matrix from  I-RF to PF-RF
rot_to_ECI = rot_to_PF';     % Coord. Transf. Matrix from PF-RF to  I-RF

% Multiply by Rotation matrix. 
r_old     = [x_orb; y_orb; z_orb]; % (!should of been a series of column vectors.)
r_new     = rot_to_ECI*r_old; 
r_current = r_new(:,end); % Last entry was based on True anomaly. coe(6)
r_new(:,end) = []; %Delete the last entry

%% Plotting

% Plot Orbit "Black" - "Orbit Initial Plot"
oiplot = plot3(UIaxishand, r_new(1,:), r_new(2,:), r_new(3,:)); % 
oiplot.Color = "k";
hold(UIaxishand, 'on'); % 

% Plot Current Position "Green" - "Orbit Initial Location Plot"
oilplot = plot3(UIaxishand, r_current(1,1), r_current(2,1), r_current(3,1)); 
oilplot.Marker          = 'o';
oilplot.MarkerEdgeColor = 'g';
oilplot.MarkerSize      = 10;

oiplot2 = [];
oilplot2 = [];
if nargin > 4 % if a second figure is specified.
% Plot Orbit "Black" - "Orbit Initial Plot"
oiplot2 = plot3(UIaxishand2, r_new(1,:), r_new(2,:), r_new(3,:)); % 
oiplot2.Color = "k";

% Plot Current Position "Green" - "Orbit Initial Location Plot"
oilplot2 = plot3(UIaxishand2, r_current(1,1), r_current(2,1), r_current(3,1)); 
oilplot2.Marker          = 'o';
oilplot2.MarkerEdgeColor = 'g';
oilplot2.MarkerSize      = 10;

else
end

fprintf('Start "%s"\n',fxname)
end
