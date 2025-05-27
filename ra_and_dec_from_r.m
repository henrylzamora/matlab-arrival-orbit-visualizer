% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [ra dec] = ra_and_dec_from_r(r,output)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{  
  This function calculates the right ascension and the
  declination from the geocentric equatorial position vector

  r       - position vector
  l, m, n - direction cosines of r
  ra      - right ascension (degrees)
  dec     - declination (degrees)
%}
% ----------------------------------------------
if nargin < 2
    output = false; % Default behavior: do not print results
end
l = r(1)/norm(r);
m = r(2)/norm(r);
n = r(3)/norm(r);

dec = asind(n);

if m > 0
    ra = acosd(l/cosd(dec));
else
    ra = 360 - acosd(l/cosd(dec));
end


if output == true % Hide printout if not desired.
fprintf('\nOutputs:"ra_and_dec_from_r.m"\n')
fprintf('     Right Ascension, ra  = %.2f degrees \n', ra)
fprintf('     Declination    , dec = %.2f degrees \n', dec)
% fprintf('\n End. \n')
end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~