function [RE,oblateness,mu,TE,J2,Color,mass] = func_Primary_Parameters(Primary)
    % Planetary Parameters
    % RE: Equatorial Radius                [km]
    % oblateness                           [-]
    % mu: Gravitational Parameter          [km^3/s^2]
    % TE: Sidereal Period                  [s]
    % J2: Second Zonal Spherical Harmonic  [-]
    % mass:                                [kg]

    
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
    else
        fprintf('Error in inputs. func_Primary_Parameters \n')
    end
end