function [ v_I_new] = func_change_velort_RevB(ra, dec, rot_axis,rot_angle, v_I)
% Function: Change Velocity Orientation.
% Function takes the geocentric Velocity, and the desired Rotation angle
% and Rotates the Velocity Vector in the Topocentric SEZ Frame. 
% Makes Rotations in the SEZ frame. 
% Transforms back to the GCI Fram.e
% Outputs the new Rotated Velocity in the Geocentric GCI Frame (IJK)
% fxname = "func_change_velort_RevB.m";  fprintf('Start "%s"\n',fxname) 

% Geocentric To Topocentric Frame Conversion
fprintf('Rotation in between GCI and SEZ Frame. Matrices\n')
    rot2 = [cosd(-dec )   0      sind(-dec );
            0             1       0;
            -sind(-dec )   0       cosd(-dec)] %3-15 vallado 5th
    
    rot3 =[ cosd(-ra)  sind(-ra)  0 
           -sind(-ra)  cosd(-ra)  0
               0       0     1]
% rot2_RevA = axang2rotm([ 0 1 0 -dec*pi/180])
% rot3_RevA = axang2rotm([ 0 0 1 -ra *pi/180])

    % SEZ    =           IJK
    v_I_sez  = rot2*rot3*v_I; % Velocity Vector
    % v_I_GCI  = rot3'*rot2'*v_I_sez;



% Rotation within SEZ Fram
    switch rot_axis % Rotation Axis
        % Case Structure to Pick Axis of Rotation
        %       rot_axis = 1: Rotated about S for Elevation Angle Change
        %                = 3:               Z for Azimuth   Angle Change
        % Input angle in Degrees
        case 1 % Rotate about S vector: (Elevation Angle.
            fprintf('Rotation in SEZ Frame. Elevation\n')
            rot1_local =[ 1         0               0
                          0     cosd(rot_angle)  sind(rot_angle) 
                          0    -sind(rot_angle)  cosd(rot_angle)]
            v_I_sez_new = rot1_local*v_I_sez; % New, Rotated Velocity, SEZ 


        case 3 % Rotate about Azimuth
            fprintf('Rotation in SEZ Frame. Azimuth\n')
            rot3_local =[ cosd(rot_angle)  sind(rot_angle)  0 
                         -sind(rot_angle)  cosd(rot_angle)  0
                                  0          0              1]
            v_I_sez_new = rot3_local*v_I_sez; % New, Rotated Velocity, SEZ 
        otherwise
            fprintf('Error in func_change_velort.m\n')
            fprintf('Case not Defined.\n')
    end

    v_I_new = rot3'*rot2'*v_I_sez_new; % New Velocity Vector, GCI
                                 % The inverse process requires revers order of
                             
fprintf('Rotation Angle: %0.3f degrees, about axis %0f\n',  rot_angle,rot_axis)
fprintf('Old Velocity, v_I      = %.3f, %.3f, %.3f  km^2/s \n', v_I)   
fprintf('New Velocity, v_I      = %.3f, %.3f, %.3f  km^2/s \n', v_I_new)

fprintf('Old Velocity, v_I_SEZ  = %.3f, %.3f, %.3f  km^2/s \n', v_I_sez)
fprintf('New Velocity, v_I_SEZ  = %.3f, %.3f, %.3f  km^2/s \n', v_I_sez_new)
% fprintf('End "%s"\n',fxname) 
