function thermCond_snow = thermCond_snow_compiled2(ice_snow, water_snow, layerThick_snow)

snow_density = max(0,min(1, (ice_snow + water_snow) ./ max(1e-20, layerThick_snow)));
thermCond_snow = max(5e-2, 2.3.*snow_density.^1.88);

