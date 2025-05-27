function [tilt,rotm] = func_obliq(Primary)
    % Planetary Parameters
    % RE: Equatorial Radius                [km]
    % oblateness                           [-]
    % mu: Gravitational Parameter          [km^3/s^2]
    % TE: Sidereal Period                  [s]
    % J2: Second Zonal Spherical Harmonic  [-]
    % mass:                                [kg]

    % tilt: planetary tilt [degrees]
    % rotm: [3 x 3 matrix ] "pre multiply"

    if strcmp(Primary,'Mercury')
        tilt = 0.01; %degrees.
        axang = [1 0 0 tilt*pi/180]; % Radians
        rotm = axang2rotm(axang);
    
    elseif strcmp(Primary,'Venus')
        tilt = 2.64; %degrees.
        axang = [1 0 0 tilt*pi/180] ;% Radians
        rotm = axang2rotm(axang);;
    
    elseif strcmp(Primary,'Earth')    
        tilt = 23.44; %degrees.
        axang = [1 0 0 tilt*pi/180]; % Radians
        rotm = axang2rotm(axang);
    
    elseif strcmp(Primary,'Mars')   
        tilt = 25.19; %degrees.
        axang = [1 0 0 tilt*pi/180]; % Radians
        rotm = axang2rotm(axang);
    
    elseif strcmp(Primary,'Jupiter')
        tilt = 3.12; %degrees.
        axang = [1 0 0 tilt*pi/180]; % Radians
        rotm = axang2rotm(axang);
    
    elseif strcmp(Primary,'Saturn') 
        tilt = 26.73; %degrees.
        axang = [1 0 0 tilt*pi/180]; % Radians
        rotm = axang2rotm(axang);
    
    elseif strcmp(Primary,'Uranus') 
        tilt = 82.23; %degrees.
        axang = [1 0 0 tilt*pi/180] ;% Radians
        rotm = axang2rotm(axang);
    
    elseif strcmp(Primary,'Neptune')
        tilt = 28.33; %degrees.
        axang = [1 0 0 tilt*pi/180]; % Radians
        rotm = axang2rotm(axang);
    
    elseif strcmp(Primary,'Moon')
        tilt = 60.41; %degrees.
        axang = [1 0 0 tilt*pi/180]; % Radians
        rotm = axang2rotm(axang);
    end
end