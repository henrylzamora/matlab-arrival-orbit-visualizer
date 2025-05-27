function func_ClearGUI_values()
% function to reset or zero out values. 
% Used when changing settings, to ensure new calculations are updated 
% appropriately.
fprintf('\nfunc_ClearGUI Values. Start\n')
app.SelectPlanetDropDown.Value = "";
app.Planet = app.SelectPlanetDropDown.Value;
app.RightAscSlider.Value    = 0;  app.DeclinationSlider.Value = 0;    % Slider Values  
app.degRightAsc.Value       = 0;  app.degDeclin.Value         = 0;    % Readout Display

app.El_elSlider.Value       = 0;  app.Az_betaSlider.Value       = 0;  % Slider Values
app.degEl_el.Value          = 0;  app.degAz_beta.Value          = 0;  % Readout Display
end
