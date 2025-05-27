function [azm, elev,v_I_sez ] = func_GCI_to_SEZ_velort_RevB(dec, ra, v_I)
% Function outputs the Azimuth and Elevation angles in the Topocentric
% reference frame (SEZ: South, East, Zenith)
% Geocentric To Topocentric
            % rot1 = [1       0          0
            %         0   cosd(dec)  sind(dec)
            %         0  -sind(dec)  cosd(dec)];
% fxname = "func_GCI_to_SEZ_velort_RevB.m";
% fprintf('Start "%s"\n',fxname) 
arg = -dec;
rot2 = [cosd(arg )   0      sind(arg )
        0             1       0
        -sind(arg )   0       cosd(arg)]; %3-15 vallado 5th

rot3 =[ cosd(-ra)  sind(-ra)  0 
       -sind(-ra)  cosd(-ra)  0
           0       0     1];
% v_I
v_I_sez  = rot2*rot3*v_I;
visezmag = norm(v_I_sez);
el       = asind(v_I_sez(3)/visezmag);
el_I     = asind(v_I(3)/norm(v_I ));
if abs(el) ==90;
            % azm  = asind( v_I_sez(2)/(sqrt(v_I_sez(1)^2+sqrt(v_I_sez(2)^2));
            % azmb = asind(-v_I_sez(1)/(sqrt(v_I_sez(1)^2+sqrt(v_I_sez(2)^2));
    fprintf('Undefined case. Velocity normal to surface. \n')
    azm = NaN;
    return
else
    % SEZ
    den  = sqrt(  v_I_sez(1)^2 + v_I_sez(2)^2 );
            % azm  = asind( v_I_sez(2)/ den);
            % azmb = acosd(-v_I_sez(1)/ den);
    azm = atan2d(v_I_sez(2),-v_I_sez(1));

    % % Geo Centric Inertial
    % den2   = sqrt( v_I(1)^2 + v_I(2)^2 );
    % azm_I  = asind( v_I(2)/ den2);
    % azmb_I = acosd(-v_I(1)/ den2);
end

elev = el;

fprintf('     Elevation and Azimuth of V_inf Vector\n')
% fprintf('     IJK: Elevation angle: %0.2f degrees\n', el_I)
% fprintf('     IJK: Azimuth angle:  asin() =beta = %0.2f degrees\n', azm_I)

fprintf('     SEZ: Elevation angle: %0.2f degrees\n', el)
fprintf('     SEZ: Azimuth angle:  asin() =beta = %0.2f degrees\n\n', azm)
end


