
function [surf_T, melt, sublimation, snowfall, rainfall, T, d_energy, d_new_snow_layerThick, d_new_melt_layerThick] = get_boundary_condition_u_compiled2(T, surf_T, melt, sublimation,...
 day_sec, timestep, snowfall, rainfall, ice_snow, upper_cell, d_energy, snow_mat1, L_f, c_i, layerThick_snow, wind_speed_class)


%in: tile.FORCING.TEMP.surfT, tile.ENSEMBLE.PARA.ensemble_size,
%tile.FORCING.TEMP.melt_bare, tile.ENSEMBLE.STATVAR.melt_fraction, 
%tile.FORCING.TEMP.melt_forest, ground.CONST.day_sec, tile.timestep,
%tile.FORCING.TEMP.snowfall, tile.ENSEMBLE.STATVAR.snowfall_factor, 
%ground.STATVAR.ice_snow, ground.STATVAR.upper_cell,
%ground.TEMP.d_energy,ground.TEMP.snow_mat1, ground.CONST.L_f,
%ground.CONST.c_i, ground.STATVAR.layerThick_snow, tile.ENSEMBLE.STATVAR.wind_speed_class

%out: : surf_T, melt, snowfall, T, TEMP: d_energy, d_new_snow_layerThick


                melt = max(0, melt ./ 1000 ./day_sec .* timestep);  %in [m], constant timestep
                sublimation = sublimation ./ 1000 ./ day_sec .* timestep;
                snowfall = snowfall ./1000 ./ day_sec .* timestep;
                rainfall = rainfall ./1000 ./ day_sec .* timestep;

                melt = max(0, min(melt, sum(ice_snow,1) + snowfall - sublimation)); %limit the melt, so that it doesn't exceed the existing snow
                surf_T = double((snowfall  - sublimation - melt + sum(ice_snow, 1)) <=0 | (surf_T <0)) .* surf_T; %set to zero if there is snow and T is positive
                surf_T = double(rainfall < 5) .* surf_T;
                snowfall = snowfall + rainfall .* double(rainfall < 5);
                rainfall = rainfall - rainfall .* double(rainfall < 5);
                for i=1:4 %assign boundary condition T to correct cell and all cells above
                    T(i, :) =  T(i, :) + double(i<=upper_cell) .* (surf_T - T(i, :));
                end

                d_energy(1:4,:) = d_energy(1:4,:) + snow_mat1 .*(repmat(snowfall./timestep .* ...
                    (-L_f + surf_T .* c_i), 4, 1)  - repmat((melt + sublimation)./ timestep, 4, 1) ...
                    .* (-L_f + T(2:5,:) .* c_i));

                a = 109;
                b = 6;
                c = 26;
                Ts=min(0,surf_T);
                new_snow_density = max(50, a + b.*Ts + (c .* wind_speed_class.^0.5));
                %melt > snowfall: no increase with new snow density (all new snow melts), decrease with existing snow density
                %melt < snowfall: no net melt,

                d_new_snow_layerThick = repmat(double(melt + sublimation < snowfall) .* (snowfall - melt - sublimation) .* 920 ./ new_snow_density , 4, 1) .* snow_mat1;
                d_new_melt_layerThick = repmat(double(melt + sublimation > snowfall) .* (melt + sublimation - snowfall ), 4, 1) .* snow_mat1 .* layerThick_snow ./ max(1e-10, ice_snow);
