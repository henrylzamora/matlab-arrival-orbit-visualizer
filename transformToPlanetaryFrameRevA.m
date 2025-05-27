function [ rotma] = transformToPlanetaryFrameRevA( alpha_deg, delta_deg)
    % Function to transform coordinates from ICRF to planetary body-inertial frame
    % Input:
    %   vec_old    - 3x1 vector in initial coordiate reference system
    %   alpha_deg  - Right ascension of planet's north pole (degrees)
    %   delta_deg  - Declination of planet's north pole (degrees)

    % Output:
    %   vec_new     - 3x1 vector in planetary body-fixed coordinates

    % Convert angles to radians
    alpha = deg2rad(alpha_deg);
    delta = deg2rad(delta_deg);

    % % Rotation matrices
    % R1 = [cos(W)  sin(W) 0;
    %      -sin(W)  cos(W) 0;
    %            0      0  1];

    R2 = [1     0          0;
          0 cos(delta) sin(delta);
          0 -sin(delta) cos(delta)];

    R3 = [cos(alpha) sin(alpha) 0;
         -sin(alpha) cos(alpha) 0;
               0          0     1];

    % Full rotation matrix: from ICRF to planet-fixed frame
    rotma = R2 * R3;

    % % Apply transformation
    % vec_new = rotma * vec_old;
end