function [r_entry] = Update_R_Entry (ROI,declin,ra)
% Declination is theta
% Right Ascencion is phi.
try
xp  = ROI*cosd(declin)*cosd(ra);     % position components, xp = xpoint
yp  = ROI*cosd(declin)*sind(ra); 
zp  = ROI*sind(declin);
r_entry = [xp,yp,zp];
catch
whos ROI declin ra
ROI
end

end