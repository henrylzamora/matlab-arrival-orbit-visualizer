function [alpha0, delta0, W] = getRotationalElements(planet, JD)
    % Function to return rotational elements (alpha, delta, W) for planets
    % Input:
    %   planet - Name of the planet as a string (e.g., 'Mercury', 'Mars')
    %   JD - Julian date for computation
    % Output:
    %   alpha - Right Ascension of north pole (degrees)
    %   delta - Declination of north pole (degrees)
    %   W - Prime meridian angle (degrees)

    % Define constants
    J2000 = 2451545.0; % Standard epoch
    d = JD - J2000;    % Days since J2000
    T = d / 36525;     % Julian centuries since J2000

 if strcmp(planet,'Sun')
    alpha0 = 286.13;
    delta0 = 63.87;
    W = 84.176 + 14.1844000*d; %(a)
elseif strcmp(planet,'Mercury')
    
    alpha0 = 281.0103 - 0.0328 *T;
    delta0 = 61.4155 - 0.0049 *T;
    % Define
    M1 = 174.7910857 + 4.092335*d;
    M2 = 349.5821714 + 8.184670*d;
    M3 = 164.3732571 + 12.277005*d;
    M4 = 339.1643429 + 16.369340*d;
    M5 = 153.9554286 + 20.461675*d;
    % Assume Plus↓
    W = 329.5988 + 0.0037 + 6.1385108*d...
    +0.01067257 *sin(M1)...
    -0.00112309 *sin(M2)...
    -0.00011040 *sin(M3)...
    -0.00002539 *sin(M4)...
    -0.00000571 *sin(M5);...


elseif strcmp(planet,'Venus')
    alpha0 = 272.76;%(b)
    delta0 = 67.16;
    W = 160.20 - 1.4813688*d;%(c)


elseif strcmp(planet,'Mars') 
alpha0 = 317.269202 - 0.10927547*T...
    +0.000068 *sin(198.991226 + 19139.4819985*T )...
    +0.000238 *sin(226.292679 + 38280.8511281*T )...
    +0.000052 *sin(249.663391 + 57420.7251593*T )...
    +0.000009 *sin(266.183510 + 76560.6367950*T )...
    +0.419057 *sin(79.398797 + 0.5042615*T );
delta0 = 54.432516 - 0.05827105*T...
    +0.000051 *cos(122.433576 + 19139.9407476*T )...
    +0.000141 *cos(43.058401 + 38280.8753272*T )...
    +0.000031 *cos(57.663379 + 57420.7517205*T )...
    +0.000005 *cos(79.476401 + 76560.6495004*T )...
    +1.591274 *cos(166.325722 + 0.5042615*T );
W = 176.049863 + 350.891982443297*d...
    +0.000145 *sin(129.071773 + 19140.0328244*T )...
    +0.000157 *sin(36.352167 + 38281.0473591*T )...
    +0.000040 *sin(56.668646 + 57420.9295360*T )...
    +0.000001 *sin(67.364003 + 76560.2552215*T )...
    +0.000001 *sin(104.792680 + 95700.4387578*T )...
    +0.584542 *sin(95.391654 + 0.5042615*T );%(d)
elseif strcmp(planet,'Jupiter') 
    
    Ja = 99.360714 + 4850.4046*T; Jb = 175.895369 + 1191.9605*T;
    Jc = 300.323162 + 262.5475*T; Jd = 114.012305 + 6070.2476*T;
    Je = 49.511251 + 64.3000*T;
    
    alpha0 = 268.056595 - 0.006499*T + 0.000117 *sin(Ja) + 0.000938 *sin(Jb)...
        +0.001432 *sin(Jc) + 0.000030 *sin (Jd) + 0.002150 *sin (Je);
    delta0 = 64.495303 + 0.002413*T + 0.000050 *cos(Ja) + 0.000404*cos(Jb)...
        +0.000617 *cos(Jc) - 0.000013 *cos(Jd) + 0.000926 *cos(Je);
    W = 284.95 + 870.5360000*d;%(e)

elseif strcmp(planet,'Saturn')  
    alpha0 = 40.589 - 0.036*T;
    delta0 = 83.537 - 0.004*T;
    W = 38.90 + 810.7939024*d; %(e)
elseif strcmp(planet,'Uranus')   
   alpha0 = 257.311;
    delta0 = -15.175;
    W = 203.81 - 501.1600928*d;%(e);
elseif strcmp(planet,'Neptune') 
     N = 357.85 + 52.316*T %(f)
    alpha0 = 299.36 + 0.70 *sin(N);
    delta0 = 43.46 - 0.51 *cos(N);
    W = 249.978 + 541.1397757*d - 0.48 *sin(N);
   
 else
     fprintf('Celestial Body not supported \n')
     return
 end
 W = mod(W,360); 
% Notes
% (a) The equation W for the Sun is now corrected for light travel time and removing the aberration correction.
%     See the Appendix in Seidelmann et al. (2007)
% (b) The 20◦ meridian of Mercury is defined by the crater Hun Kal
% (c) The 0◦ meridian of Venus is defined by the central peak in the crater Ariadne
% (d) The longitude of the Viking 1 lander on Mars is defined to be 47◦.95137 west (Kuchynka et al. 2014), maintaining the 0◦ meridian through the crater Airy-0
% (e) The equations for W for Jupiter, Saturn, and Uranus refer to the rotation of their magnetic fields (System
%      III). On Jupiter, System I (WI = 67◦.1 + 877◦.900d) refers to the mean atmospheric equatorial rotation;
%     System II (WI I = 43◦.3+870◦.270d) refers to the mean atmospheric rotation north of the south component of the north equatorial belt, and south of the north component of the south equatorial belt
% (f) The equations for Neptune refer to the rotation of optically observed features in the Neptunian atmosphere
%     (System II), while still using the previous expressions for pole position and precession

end
