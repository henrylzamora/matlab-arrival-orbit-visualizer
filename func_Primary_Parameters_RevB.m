function varargout = func_Primary_Parameters(Primary)
% function [RE,oblateness,mu,TE,J2,Color,mass] = func_Primary_Parameters(Primary)
% function [RE,RP,mu,TE,J2,Color_Primary,Color_Primary_Neg,Map_Primary,r_SOI] = func_Primary_Parameters1(Primary)
% function [planetparam]

if nargout==1 % Planet parameters as a structure
    if strcmp(Primary,'Mercury')
        planetparam.RE = 2440;
        planetparam.mu = 22032;
        planetparam.TE = 58.65*86400;
        planetparam.J2 = 60E-6;
        planetparam.oblateness = 0;
        planetparam.Color = [.55,.55,.56];
        planetparam.mass = 330.2 * 10^21;
        f = 0;
        r = 57.91e6;
    elseif strcmp(Primary,'Venus')
        planetparam.RE = 6052;
        planetparam.mu = 324859;
        planetparam.TE = -243*86400;
        planetparam.J2 = 4.458E-6;
        planetparam.oblateness = 0;
        planetparam.Color = [1,.5,0];
        planetparam.mass = 4.869 * 10^24;
        f = 0;
        r = 108.2e6;

    elseif strcmp(Primary,'Earth')    
        planetparam.RE = 6378;
        planetparam.mu = 398600;
        planetparam.TE = 23.9345*3600;
        planetparam.J2 = 1.08263E-3;
        planetparam.oblateness = 0.003353;
        planetparam.Color = [.2,.2,1];
        planetparam.mass = 5.974 * 1024;
        f = 0.0033528131;
        r = 149.6E6;
    elseif strcmp(Primary,'Mars')   
        planetparam.RE = 3396;
        planetparam.mu = 42828;
        planetparam.TE = 24.62*3600;
        planetparam.J2 = 1.96045E-3;
        planetparam.oblateness = 0.00648;
        planetparam.Color = [1,.4,.1];
        planetparam.mass = 641.9 * 10^21;
                f = 0.00648;
        Color_Primary = [1,.4,.1];
        r = 227.9e6;
    elseif strcmp(Primary,'Jupiter')
        planetparam.RE = 71490;
        planetparam.mu = 126686534;
        planetparam.TE = 9.925*3600;
        planetparam.J2 = 14.736E-3;
        planetparam.oblateness = 0.06487;
        planetparam.Color = [.9,.7,.2];
        planetparam.mass = 1.899 * 10^27
                f = 0.06487;
        Color_Primary = [.9,.7,.2];
        r = 778.6e6;
    elseif strcmp(Primary,'Saturn') 
        planetparam.RE = 60270;
        planetparam.mu = 37931187;
        planetparam.TE = 10.66*3600;
        planetparam.J2 = 16.298E-3;
        planetparam.oblateness = 0.09796;
        planetparam.Color = [1,.9,.6];
        planetparam.mass = 568.5 * 10^24;
        f = 0.09796;
        Color_Primary = [1,.9,.6];
        r = 1433e6;
    elseif strcmp(Primary,'Uranus') 
        planetparam.RE = 25560;
        planetparam.mu = 5793939;
        planetparam.TE = -17.24*3600;
        planetparam.J2 = 3.34343E-3;
        planetparam.oblateness = 0.02293;
        planetparam.Color = [0,.9,1]; 
        planetparam.mass = 86.83 * 10^24;
                f = 0.02293;
        Color_Primary = [0,.9,1]; 
        r = 2872e6;
    elseif strcmp(Primary,'Neptune')
        planetparam.RE = 24764;
        planetparam.mu = 6836529;
        planetparam.TE = 16.21*3600;
        planetparam.J2 = 3.411E-3;
        planetparam.oblateness = 0.01708;
        planetparam.Color = [0,.5,1];
        planetparam.mass = 102.4 * 10^24
                f = 0.01708;
        Color_Primary = [0,.5,1];
        r = 4495e6; 
    elseif strcmp(Primary,'Moon')
        planetparam.RE = 1737;
        planetparam.mu = 4905;
        planetparam.TE = 27.32*86400;
        planetparam.J2 = 202.7E-6;
        planetparam.oblateness = 0.0012;
        planetparam.Color = [.6,.61,.6]; 
        planetparam.mass = 73.48 * 10^21;
    end
            % Polar Radius
        planetparam.RP = planetparam.RE*(1-planetparam.f);
    
        % Colors
        HSV_CP = rgb2hsv(Color_Primary);
        if HSV_CP(2)<0.5
            if HSV_CP(1)<0.5
                Color_Primary_Neg = hsv2rgb([HSV_CP(1)+0.5,1,1]);
            else
                Color_Primary_Neg = hsv2rgb([HSV_CP(1)-0.5,1,1]);
            end
        else
            Color_Primary_Neg = [1,1,1]-Color_Primary;
        end
    
        % Planet's map
        image_file  = sprintf('Map_%s.jpg',Primary);% Changed the filename Convenction
        Map_Primary = imread(image_file);           % now "Map_planet"
    
        % Radius of Sphere of Influence (SOI)
        mu_Sun = 1.3271244e11; %[km^3/s^2]
        r_SOI = r*(mu/mu_Sun)^(2/5); %[km]

    % Planetary Parameters
    % RE: Equatorial Radius                [km]
    % oblateness                           [-]
    % mu: Gravitational Parameter          [km^3/s^2]
    % TE: Sidereal Period                  [s]
    % J2: Second Zonal Spherical Harmonic  [-]
    % mass:                                [kg]
elseif nargout == 7

    if strcmp(Primary,'Mercury')
        RE = 2440;
        mu = 22032;
        TE = 58.65*86400;
        J2 = 60E-6;
        oblateness = 0;
        Color = [.55,.55,.56];
        mass = 330.2 * 10^21;
    elseif strcmp(Primary,'Venus')
        RE = 6052;
        mu = 324859;
        TE = -243*86400;
        J2 = 4.458E-6;
        oblateness = 0;
        Color = [1,.5,0];
        mass = 4.869 * 10^24;
    elseif strcmp(Primary,'Earth')    
        RE = 6378;
        mu = 398600;
        TE = 23.9345*3600;
        J2 = 1.08263E-3;
        oblateness = 0.003353;
        Color = [.2,.2,1];
        mass = 5.974 * 1024;
    elseif strcmp(Primary,'Mars')   
        RE = 3396;
        mu = 42828;
        TE = 24.62*3600;
        J2 = 1.96045E-3;
        oblateness = 0.00648;
        Color = [1,.4,.1];
        mass = 641.9 * 10^21;
    elseif strcmp(Primary,'Jupiter')
        RE = 71490;
        mu = 126686534;
        TE = 9.925*3600;
        J2 = 14.736E-3;
        oblateness = 0.06487;
        Color = [.9,.7,.2];
        mass = 1.899 * 10^27
    elseif strcmp(Primary,'Saturn') 
        RE = 60270;
        mu = 37931187;
        TE = 10.66*3600;
        J2 = 16.298E-3;
        oblateness = 0.09796;
        Color = [1,.9,.6];
        mass = 568.5 * 10^24;
    elseif strcmp(Primary,'Uranus') 
        RE = 25560;
        mu = 5793939;
        TE = -17.24*3600;
        J2 = 3.34343E-3;
        oblateness = 0.02293;
        Color = [0,.9,1]; 
        mass = 86.83 * 10^24;
    elseif strcmp(Primary,'Neptune')
        RE = 24764;
        mu = 6836529;
        TE = 16.21*3600;
        J2 = 3.411E-3;
        oblateness = 0.01708;
        Color = [0,.5,1];
        mass = 102.4 * 10^24
    elseif strcmp(Primary,'Moon')
        RE = 1737;
        mu = 4905;
        TE = 27.32*86400;
        J2 = 202.7E-6;
        oblateness = 0.0012;
        Color = [.6,.61,.6]; 
        mass = 73.48 * 10^21;
    end
        varargout = [RE,oblateness,mu,TE,J2,Color,mass] ;
    
    
    
    elseif nargout == 9
% function [RE,RP,mu,TE,J2,Color_Primary,Color_Primary_Neg,Map_Primary,r_SOI] = func_Primary_Parameters1(Primary)
    % Planetary Parameters
    % RE: Equatorial Radius                [km]
    % RP: Polar Radius                     [km]
    % f:  flatness                         [-]
    % mu: Gravitational Parameter          [km^3/s^2]
    % TE: Sidereal Period                  [s]
    % J2: Second Zonal Spherical Harmonic  [-]
    if strcmp(Primary,'Mercury')
        RE = 2440;
        mu = 22032;
        TE = 58.65*86400;
        J2 = 60E-6;
        f  = 0;
        Color_Primary = [.55,.55,.56];
        r = 57.91e6;
    elseif strcmp(Primary,'Venus')
        RE = 6052;
        mu = 324859;
        TE = -243*86400;
        J2 = 4.458E-6;
        f = 0;
        Color_Primary = [1,.5,0];
        r = 108.2e6;
    elseif strcmp(Primary,'Earth')    
        RE = 6378.1363;
        mu = 398600.4415;
        TE = 0.99726968*86400;
        J2 = 1.0826269E-3;
        f = 0.0033528131;
        Color_Primary = [.2,.2,1];
        r = 149.6E6;
    elseif strcmp(Primary,'Mars')   
        RE = 3396;
        mu = 42828;
        TE = 24.62*3600;
        J2 = 1.96045E-3;
        f = 0.00648;
        Color_Primary = [1,.4,.1];
        r = 227.9e6;
    elseif strcmp(Primary,'Jupiter')
        RE = 71490;
        mu = 126686534;
        TE = 9.925*3600;
        J2 = 14.736E-3;
        f = 0.06487;
        Color_Primary = [.9,.7,.2];
        r = 778.6e6;
    elseif strcmp(Primary,'Saturn') 
        RE = 60270;
        mu = 37931187;
        TE = 10.66*3600;
        J2 = 16.298E-3;
        f = 0.09796;
        Color_Primary = [1,.9,.6];
        r = 1433e6;
    elseif strcmp(Primary,'Uranus') 
        RE = 25560;
        mu = 5793939;
        TE = -17.24*3600;
        J2 = 3.34343E-3;
        f = 0.02293;
        Color_Primary = [0,.9,1]; 
        r = 2872e6;
    elseif strcmp(Primary,'Neptune')
        RE = 24764;
        mu = 6836529;
        TE = 16.21*3600;
        J2 = 3.411E-3;
        f = 0.01708;
        Color_Primary = [0,.5,1];
        r = 4495e6; 
    elseif strcmp(Primary,'Moon')
        RE = 1737;
        mu = 4905;
        TE = 27.32*86400;
        J2 = 202.7E-6;
        f = 0.0012;
        Color_Primary = [.6,.61,.6];
        r = 384.4e3; 
    end
        % Polar Radius
        RP = RE*(1-f);
    
        % Colors
        HSV_CP = rgb2hsv(Color_Primary);
        if HSV_CP(2)<0.5
            if HSV_CP(1)<0.5
                Color_Primary_Neg = hsv2rgb([HSV_CP(1)+0.5,1,1]);
            else
                Color_Primary_Neg = hsv2rgb([HSV_CP(1)-0.5,1,1]);
            end
        else
            Color_Primary_Neg = [1,1,1]-Color_Primary;
        end
    
        % Planet's map
        image_file  = sprintf('Map_%s.jpg',Primary);% Changed the filename Convenction
        Map_Primary = imread(image_file);           % now "Map_planet"
    
        % Radius of Sphere of Influence (SOI)
        mu_Sun = 1.3271244e11; %[km^3/s^2]
        r_SOI = r*(mu/mu_Sun)^(2/5); %[km]
                whos RE RP mu TE J2 Color_Primary Color_Primary_Neg Map_Primary r_SOI
                  % [RE,RP,mu,TE,J2,Color_Primary,Color_Primary_Neg,Map_Primary,r_SOI]
        % varargout = {RE,RP,mu,TE,J2,Color_Primary(1),Color_Primary_Neg(1),Map_Primary(1),r_SOI};
        varargout = {RE,RP,mu,TE,J2,Color_Primary,Color_Primary_Neg,Map_Primary,r_SOI};
    else
        fprintf('Invalid number of output arguments.\n')
        fprintf('Use 7 or 9 outputs as an arrary.\n');
        fprintf('Or use 1 output as structure.')
        return
    end
    
end