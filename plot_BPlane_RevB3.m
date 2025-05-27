
function [bplot, p1,p2,p3,p4,ECI,Br,Bt,Bmag,vinf] = plot_BPlane_RevB3(r, v, Planet,UIaxhand,UIaxhand2 )
%% B Plane Parameters from R and V
% Takes a user input R and V for a Spacecraft approaching a planets Sphere of Influence (SOI). %
% Calculates and Plots Br and Bt. 
    % UIaxishand - Figure 1: Far View
    % UIaxishand - Figure 2: Near View.

% Administrative
k = 3; % Scaling factor. Used in conjunction with RE. 
fxname = "plot_BPlane_RevB3.m"; fprintf('Start "%s"\n',fxname)
s = size(r);
if s(1)==1 % If row vector
    r = r';% Make it a column vector
end

planetparam = func_Primary_Parameters_RevC(Planet);
RE = planetparam.RE;
mu = planetparam.mu;

    h       = cross(r,v);
    hhat    = cross(r,v) / norm(cross(r,v) );  %
    hnorm   = norm(cross(r,v));                 % Checked good. Matches original COE
    e       =( (norm(v)^2 - mu/norm(r))*r - dot(r,v)*v )/mu;  % Vallado 2.78  Calculates e differently 
    enorm   = norm(e); % Checked good. Matches original COE
    a       = hnorm^2 /(mu*(enorm^2 - 1));
    vinf    = sqrt(mu/a);
    
    b = hnorm^2 /(mu*sqrt(enorm^2 - 1));	% Semiminor axis.
    % b3       = a*sqrt(norm(e)^2 - 1); % Aiming radius
    phi_s = acos(1/norm(e));          % Supplementary Angle, Rads
    phi_s_deg = phi_s*180/pi;         % Supplementary Angle, Deg
    K = [0 0 1];
% B Plane Unit Vectors
    % Vectors R,S,T form the coordinate system of the B Plane
    % S is parallel to the hyperbolic velocity.
    % R and T form the plane that is pierced by the spacecraft
    S_hat = (e/norm(e)*cos(phi_s)) + (cross(hhat, e)/norm(cross(hhat,e)))*sin(phi_s);
    T_hat = cross(S_hat,K)/norm(cross(S_hat,K));
    R_hat = cross(S_hat,T_hat);

    B     = b*cross(S_hat,hhat);     
    Bt    = dot(B,T_hat);     
    Br    = dot(B,R_hat);
    Bmag = norm(B);
    
% Function Outputs
fprintf('Vallado. B Plane Parameters \n')
fprintf('     Angular Momentum Vector, h   = %.2f, %.2f, %.2f  km^2/s \n', h)
fprintf('     Angular Momentum       , h   = %.2f km^2/s \n', hnorm)
fprintf('     Velocity Vector,         v   = %.2f, %.2f, %.2f  km^2/s \n', v)
fprintf('     Eccentricity Vector,     e   = %.2f, %.2f, %.2f, \n'     ,e )
fprintf('     Supplimentary Angle,   phi_s =  %.2f rad = %.1f degs\n'  , phi_s, phi_s_deg )
fprintf('     Eccentricity,     e   = %.2f \n'       ,enorm)
fprintf('     Aiming Radius,    b   = %.0f km \n'    ,b) 
% fprintf('                  ,    b3  = %.0f km \n'    ,b3) 
fprintf('     Magnitude of B,   B   = %.0f km \n'  ,norm(B))
fprintf('     B_radial,         B_r = %.0f km \n',Br)
fprintf('     B_tangential,     B_t = %.0f km \n',Bt)


fprintf('     Coordinate Vectors, R   = %.2f, %.2f, %.2f  \n'    , R_hat')
fprintf('     Coordinate Vectors, S   = %.2f, %.2f, %.2f   \n'   , S_hat')
fprintf('     Coordinate Vectors, T   = %.2f, %.2f, %.2f   \n'   , T_hat')
% a3      = mu/vinf^2         % Hyperbolic Excess Velocity




%% Plotting
% Plot the Equatorial Plane
    X =linspace(-RE*k,RE*k,51);
    Y =linspace(-RE*k,RE*k,51);
    [XECI,YECI]=meshgrid(X,Y);
    ZECI=0*XECI+0*YECI;
    ECI  = surf(UIaxhand , XECI,YECI,ZECI); hold(UIaxhand , 'on'); 
    ECI2 = surf(UIaxhand2, XECI,YECI,ZECI); hold(UIaxhand2, 'on'); % 
        % Formatting
        ECI.FaceColor = [1,0.2,0.2]; ECI2.FaceColor = [1,0.2,0.2];%Red!
        ECI.EdgeColor = 'none';      ECI2.EdgeColor = 'none';
        ECI.FaceAlpha = 0.3;         ECI2.FaceAlpha = 0.3;
        hold(UIaxhand, 'on');      hold(UIaxhand2, 'on'); %need only to invoke once per axis handle 


% Plot the B Plane.
    len = linspace(-RE*k, RE*k, 50);
    normal = cross(R_hat, T_hat);       % Define the normal vector to the plane
    [x, y] = meshgrid(len,len); % Generate grid points for plotting
    % Compute z values for the plane using the equation of the plane
    z = (-normal(1)*x - normal(2)*y) / normal(3);  

    bplot = surf(UIaxhand,  x, y, z);      % Plot the plane
    bplot2= surf(UIaxhand2, x, y, z);      % Plot the plane

        bplot.FaceColor = [0,0.5,1];bplot2.FaceColor = [0,0.5,1];%Blue
        bplot.EdgeColor = 'none';   bplot2.EdgeColor = 'none';
        bplot.FaceAlpha = 0.3;      bplot2.FaceAlpha = 0.3;



% Plot the vectors R, S, T, B. Zoomed Out Plot
    lim = RE*k; lw = 1.5;
    p1 = line(UIaxhand2, [0 R_hat(1)*lim], [0 R_hat(2)*lim], [0 R_hat(3)*lim],'linewidth',lw, 'Color','r'); 
    p3 = line(UIaxhand2, [0 T_hat(1)*lim], [0 T_hat(2)*lim], [0 T_hat(3)*lim],'linewidth',lw, 'Color','r'); 
     try
    p2 = line(UIaxhand2, [0 S_hat(1)*lim], [0 S_hat(2)*lim], [0 S_hat(3)*lim],'linewidth',lw, 'Color','r'); 
    p4 = line(UIaxhand2, [0 B(1)*lim], [0 B(2)*lim], [0 B(3)*lim],'linewidth',lw, 'Color','b');
    catch
    fprintf('Plotting B Plane Vectors Failed. \n')
     end
     lu1 = R_hat*lim*0.9;
     lu2 = S_hat*lim*0.9;
     lu3 = T_hat*lim*0.9;
     lu4 = B*0.95;
     
     % Annotate vectors with variable. 
text(UIaxhand2,lu1(1), lu1(2), lu1(3)*1.05,'$\mathbf{\hat{R}}$','Color','r','FontSize',14,'interpreter','latex');
text(UIaxhand2,lu2(1), lu2(2), lu2(3)*1.05,'$\mathbf{\hat{S}}$','Color','r','FontSize',14,'interpreter','latex');
text(UIaxhand2,lu3(1), lu3(2), lu3(3)*1.15,'$\mathbf{\hat{T}}$','Color','r','FontSize',14,'interpreter','latex');
text(UIaxhand2,lu4(1), lu4(2), lu4(3)*1.15,'$\mathbf{\hat{B}}$','Color','b','FontSize',14,'interpreter','latex');
fprintf('End "%s"\n',fxname)
%{
Executes Algorithm 78 (Vallado)
D. A. Vallado. Fundamentals of Astrodynamics and Applications. Microcosm
Press, fourth edition, 2013.

Takes a user input R and V and V_inf for a Spacecraft approaching Mars. 
Calculates Br and Bt
phi_s = supplementary turn angle. 


%}


end







