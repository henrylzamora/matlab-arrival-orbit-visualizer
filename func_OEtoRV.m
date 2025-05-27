function [r_I,v_I,r_PF,v_PF] = func_OEtoRV(e,a,i,argp,RAAN,TA,mu)

% TA, RAAN, i, argp - deg


    % From Orbital Elements to Position (R) and Velocity (V)
    % in ECI (_I) and Perifocal (_PF) coordinates

    norm_h = sqrt(a*mu*(e^2-1));      % Orbital Angular momentum (magnitude) (2.71% updated 10/3
    
    norm_r = norm_h^2/mu/(1+e*cosd(TA));      % Position Radius (magnitude)
    
    r_PF = norm_r*[cosd(TA),sind(TA),0]';       % Position Vector in PF-RF coords.
    v_PF = mu/norm_h*[-sind(TA),e+cosd(TA),0]'; % Velocity Vector in PF-RF coords.
    
    RO = [cosd(RAAN) sind(RAAN) 0; -sind(RAAN) cosd(RAAN) 0; 0 0 1];
    Ri = [1,0,0; 0 cosd(i) sind(i); 0 -sind(i) cosd(i)];
    Rw = [cosd(argp) sind(argp) 0; -sind(argp) cosd(argp) 0; 0 0 1];
    
    R_PI = Rw*Ri*RO;  % Coord. Transf. Matrix from  I-RF to PF-RF
    R_IP = R_PI';     % Coord. Transf. Matrix from PF-RF to  I-RF
    
    r_I = R_IP*r_PF; % Position Vector in I-RF coords.
    v_I = R_IP*v_PF; % Velocity Vector in I-RF coords.

    % whos r_I v_I

end

