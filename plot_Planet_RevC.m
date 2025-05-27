function [pplot, splot, eplot, r_entry] = plot_Planet_RevC(Planet,dec,ra,UIaxishand,UIaxishand2)
% Theta and phi are in degrees.
% This function will plot a planet, a sphere of Influence, and a point on 
% the spheres surface to show s/c entry point
% UIaxishand - Figure 1: Far View
% UIaxishand - Figure 2: Near View.

% Planet data (Radius of influence and RE)
planetparam= func_Primary_Parameters_RevC(Planet);
ROI = planetparam.r_SOI;
RE  = planetparam.RE;
fprintf('Start "plot_Planet_RevC.m"\n' )
fprintf('ROI of %s, = %.0f  kms \n',Planet, ROI)

% Entry position components, xe = x entry
xe  = ROI*cosd(dec)*cosd(ra);     
ye  = ROI*cosd(dec)*sind(ra); 
ze  = ROI*sind(dec);
r_entry = [xe,ye,ze];

% Sphere Data
[x, y, z] = sphere(24); 
% Scale the sphere by the radius
xs = x * ROI;    ys = y * ROI;    zs = z * ROI;

% Plot the Sphere of Influence.    
splot = surf(UIaxishand, xs, ys, zs);
    splot.FaceColor = [0.62,0.62,0.62]; % Gray
    splot.EdgeColor = 'none';
    splot.FaceAlpha = 0.05;
    hold(UIaxishand, 'on'); % Need only to invoke once per axis handle.

    % Add latitude and longitude lines
    numLines = 24;
    dec = linspace(0, 2*pi, numLines); % Longitude angles
    ra = linspace(-pi/2, pi/2, numLines); % Latitude angles
    Color_LiGray = [0.9, 0.9, 0.9]; % Light gray
    
    % Longitude lines
    for k = 1:numLines
        xLine = ROI* cos(dec(k)) * cos(ra);
        yLine = ROI* sin(dec(k)) * cos(ra);
        zLine = ROI* sin(ra);
        plot3(UIaxishand, xLine, yLine, zLine, 'Color', Color_LiGray, 'LineStyle', '-'); % Dashed black lines for longitude
    end
    % Latitude lines
    for k = 1:numLines
        xLine = ROI* cos(dec) * cos(ra(k));
        yLine = ROI* sin(dec) * cos(ra(k));
        zLine = ROI* sin(ra(k)) * ones(size(dec));
        plot3(UIaxishand, xLine, yLine, zLine, 'Color', Color_LiGray, 'LineStyle', '-'); % Dashed black lines for latitude
    end


% Plot the Planet
xp = x * RE ;    yp = y * RE ;    zp = z * RE ;
pplotdata = {xp,yp,zp};
image_file  = sprintf('Map_%s.jpg',Planet); % Format Name for Picture file
cdata = imread(image_file); % Load Pic 
pplot = surface(UIaxishand, xp, yp, zp);
    pplot.FaceColor = 'texturemap';
    pplot.FaceAlpha = 1.0;
    pplot.EdgeColor = 'none';
    pplot.CData     = cdata;

% Plot the entry point
eplot = plot3(UIaxishand,xe, ye, ze, 'Marker', 'o','MarkerEdgeColor','m','MarkerSize',9);

% Plot Unit Vectors and Labels
lw = 1.5; % linewidth setting
iplot = line(UIaxishand, [0 0.5*ROI], [0    0   ], [0 0      ],'linewidth',lw, 'Color','k'); 
jplot = line(UIaxishand, [0     0  ], [0 0.5*ROI], [0 0      ],'linewidth',lw, 'Color','k'); 
kplot = line(UIaxishand, [0     0  ], [0    0   ], [0 0.5*ROI],'linewidth',lw, 'Color','k'); 

text(UIaxishand,0.7*ROI,0, 0.05,'$\mathbf{\hat{I}}$','Color','k','FontSize',14,'interpreter','latex');
text(UIaxishand,0,0.7*ROI, 0.05,'$\mathbf{\hat{J}}$','Color','k','FontSize',14,'interpreter','latex');
text(UIaxishand,0,0.05, 0.7*ROI,'$\mathbf{\hat{K}}$','Color','k','FontSize',14,'interpreter','latex');

if nargin > 4 % if a second figure is specified.
    
% Plot the planet
hold(UIaxishand2, 'on'); 
pplot = surface(UIaxishand2, xp, yp, zp);
    pplot.FaceColor = 'texturemap';
    pplot.FaceAlpha = 1.0;
    pplot.EdgeColor = 'none';
    pplot.CData     = cdata;

% Plot Unit Vectors and Labels
iplot = line(UIaxishand2, [0 3*RE], [0 0   ], [0 0   ],'linewidth',lw, 'Color','k'); 
jplot = line(UIaxishand2, [0  0  ], [0 3*RE], [0 0   ],'linewidth',lw, 'Color','k');
kplot = line(UIaxishand2, [0  0  ], [0 0   ], [0 3*RE],'linewidth',lw, 'Color','k'); 
text(UIaxishand2,3.5*RE,0, 0.05,'$\mathbf{\hat{I}}$','Color','k','FontSize',12,'interpreter','latex');% [0, 0.44, 0.75]
text(UIaxishand2,0,3.5*RE, 0.05,'$\mathbf{\hat{J}}$','Color','k','FontSize',12,'interpreter','latex');
text(UIaxishand2,0,0.05, 3.5*RE,'$\mathbf{\hat{K}}$','Color','k','FontSize',12,'interpreter','latex');

else
end

end


