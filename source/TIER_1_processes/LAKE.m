%========================================================================
% CryoGrid TIER1 library class for functions related to water bodies 
% S. Westermann, October 2020
%========================================================================

classdef LAKE < BASE

    methods
        
        %move melted grid cells below ice cover (frozen grid cells)
        function ground = move_ice_up(ground)
            fully_melted = (ground.STATVAR.energy >= 0);
            ground.STATVAR.energy = reorganize_cells_frozen_melted(ground, ground.STATVAR.energy, fully_melted);
            ground.STATVAR.waterIce = reorganize_cells_frozen_melted(ground, ground.STATVAR.waterIce, fully_melted);
            ground.STATVAR.mineral = reorganize_cells_frozen_melted(ground, ground.STATVAR.mineral, fully_melted);
            ground.STATVAR.organic = reorganize_cells_frozen_melted(ground, ground.STATVAR.organic, fully_melted);
            %ground.STATVAR.air = reorganize_cells_frozen_melted(ground, ground.STATVAR.air, fully_melted);
            ground.STATVAR.layerThick = reorganize_cells_frozen_melted(ground, ground.STATVAR.layerThick, fully_melted);
        end
        
        function reorganized = reorganize_cells_frozen_melted(ground, variable, fully_melted)
            melted_cells = variable(fully_melted); 
            frozen_cells = variable(~fully_melted);
            reorganized = [frozen_cells; melted_cells];
        end
        
        function ground = move_ice_up2(ground) %puts all the ice in a single slab, avoiding multiple grid cells with non-zero ice contents
            fully_melted = (ground.STATVAR.energy >= 0);
            fully_frozen = ground.STATVAR.energy <= -ground.STATVAR.waterIce .* ground.CONST.L_f;
            ground.STATVAR.energy = reorganize_cells_frozen_melted2(ground, ground.STATVAR.energy, fully_melted, fully_frozen);
            ground.STATVAR.waterIce = reorganize_cells_frozen_melted2(ground, ground.STATVAR.waterIce, fully_melted,fully_frozen);
            ground.STATVAR.mineral = reorganize_cells_frozen_melted2(ground, ground.STATVAR.mineral, fully_melted, fully_frozen);
            ground.STATVAR.organic = reorganize_cells_frozen_melted2(ground, ground.STATVAR.organic, fully_melted, fully_frozen);
            %ground.STATVAR.air = reorganize_cells_frozen_melted(ground, ground.STATVAR.air, fully_melted);
            ground.STATVAR.layerThick = reorganize_cells_frozen_melted2(ground, ground.STATVAR.layerThick, fully_melted, fully_frozen);
            
            freezing = find(ground.STATVAR.energy > -ground.STATVAR.waterIce .* ground.CONST.L_f & ground.STATVAR.energy < 0);
            if ~isempty(freezing)
                energy_freezing = sum(ground.STATVAR.energy(freezing,1));
                for i=freezing(1):freezing(end)
                    ground.STATVAR.energy(i,1) = max(-ground.STATVAR.waterIce(i,1) .* ground.CONST.L_f, energy_freezing);
                    energy_freezing = min(0,energy_freezing - ground.STATVAR.energy(i,1));
                end
            end
        end
        
        function reorganized = reorganize_cells_frozen_melted2(ground, variable, fully_melted, fully_frozen)
            melted_cells = variable(fully_melted); 
            frozen_cells = variable(fully_frozen);
            freezing_cells = variable(~fully_melted & ~fully_frozen);
            reorganized = [frozen_cells; freezing_cells; melted_cells];
        end
        
        
        %stratifies the water column, move 1/2 cell up per timestep when water density is lower
        function ground = stratify_simple(ground)  
            if size(ground.STATVAR.energy,1) > 1
                density_water = water_density(ground);
                swap = double(ground.STATVAR.energy(1:end-1,1) >=0 & ground.STATVAR.energy(2:end,1) >= 0 & density_water(2:end,1) < density_water(1:end-1,1));
                
                %only energy is moved, water not moved physically, change if there are solutes, etc.
                energy_down = swap .* 0.5 .* min(ground.STATVAR.waterIce(1:end-1,1), ground.STATVAR.waterIce(2:end,1)) ./ ground.STATVAR.waterIce(1:end-1,1) .* ground.STATVAR.energy(1:end-1,1);
                energy_up = swap .* 0.5 .* min(ground.STATVAR.waterIce(1:end-1,1), ground.STATVAR.waterIce(2:end,1)) ./ ground.STATVAR.waterIce(2:end,1) .* ground.STATVAR.energy(2:end,1);
                
                ground.STATVAR.energy(1:end-1,1) = ground.STATVAR.energy(1:end-1,1) - energy_down + energy_up;
                ground.STATVAR.energy(2:end,1) = ground.STATVAR.energy(2:end,1)  + energy_down - energy_up;
            end
        end
        
        function ground = stratify_surface_mixed_layer(ground, tile)
            density_water = water_density_UNESCO(ground, ground.STATVAR.T);
            
            top_all_water_cell = find(ground.STATVAR.energy > -ground.STATVAR.waterIce .* ground.CONST.L_f, 1, 'first'); %first cell that is not fully frozen - must be modified for saline lake
            bottom_all_water_cell = find(ground.STATVAR.energy >= 0, 1, 'last');
            if ~isempty(top_all_water_cell) && ~isempty(bottom_all_water_cell) && ~(top_all_water_cell==bottom_all_water_cell)
                
                T = ground.STATVAR.T;
                layerThick = ground.STATVAR.layerThick;
                elevations = -cumsum([0;layerThick]);
                elevations_midpoint = -(cumsum(layerThick) - layerThick./2);
                energy = ground.STATVAR.energy;
                energy_frozen = double(energy<0) .* energy; %the frozen part of the energy (unfrozen part is zero)
                layerThick_frozen = max(0, min(energy_frozen ./ (-ground.STATVAR.waterIce .* ground.CONST.L_f) .* layerThick, layerThick)); %frozen part of grid cell
                
                wind = ground.TEMP.wind .* double(energy(1,1) >= -0.02 .* ground.STATVAR.area(1,1) .* ground.CONST.L_f); %1cm ice as threshold
                energy_avail_mix = ground.TEMP.energy_avail_mix;
                
                coeff_mix_conv = ground.CONST.coeff_mix_conv;
                coeff_wind_stir = ground.CONST.coeff_wind_stir;
                coeff_mix_turb = ground.CONST.coeff_mix_turb;
                coef_wind_drag = ground.CONST.coef_wind_drag;
                g = ground.CONST.g;
                timestep = tile.timestep;

                mixed_layer_bottom_cell = top_all_water_cell;
                moment0 = 0;
                moment1 = 0;
                energy_mixing = 0;
                density_mixed_layer = density_water(mixed_layer_bottom_cell,1);
                
                while mixed_layer_bottom_cell < bottom_all_water_cell && density_mixed_layer >= density_water(mixed_layer_bottom_cell+1,1)
                    
                    density_x_lt = density_water(mixed_layer_bottom_cell,1) .* layerThick(mixed_layer_bottom_cell,1);
                    moment0 = moment0 + density_x_lt;
                    moment1 = moment1 + density_x_lt .* elevations_midpoint(mixed_layer_bottom_cell,1);
                    
                    mixed_layer_bottom_cell = mixed_layer_bottom_cell +1;
                    range = top_all_water_cell:mixed_layer_bottom_cell;
                    T_mixed_layer = sum(T(range,1) .* (layerThick(range,1)-layerThick_frozen(range,1))) ./ sum(layerThick(range,1)-layerThick_frozen(range,1));
                    density_mixed_layer = water_density_UNESCO(ground, T_mixed_layer);
                    
                end
                
                %last layer
                density_x_lt = density_water(mixed_layer_bottom_cell,1) .* layerThick(mixed_layer_bottom_cell,1);
                moment0 = moment0 + density_x_lt;
                moment1 = moment1 + density_x_lt .* elevations_midpoint(mixed_layer_bottom_cell,1);
                
                %mix the layer
                range = top_all_water_cell:mixed_layer_bottom_cell;
                T_mixed_layer = sum(T(range,1) .* (layerThick(range,1)-layerThick_frozen(range,1))) ./ sum(layerThick(range,1)-layerThick_frozen(range,1));
                T(range,1) = T_mixed_layer;
                density_mixed_layer = water_density_UNESCO(ground, T_mixed_layer);
                density_water(range,1) = density_mixed_layer;
                energy_mixed_layer = sum(energy(range,1)-energy_frozen(range,1)) ./ sum(layerThick(range,1)-layerThick_frozen(range,1)) .* (layerThick(range,1)-layerThick_frozen(range,1));
                energy(range,1) = energy_mixed_layer +  energy_frozen(range,1);
                T(range,1) = double(energy(range,1) >0) .* T(range,1);
                energy_frozen = double(energy<0) .* energy; %the frozen part of the energy (unfrozen part is zero)
                layerThick_frozen = max(0, min(energy_frozen ./ (-ground.STATVAR.waterIce .* ground.CONST.L_f) .* layerThick, layerThick)); %frozen part of grid cell
                
                if mixed_layer_bottom_cell < bottom_all_water_cell %not bottom reached
                    energy_convection = max(0, 0.5 .* g .* coeff_mix_conv ./density_mixed_layer .* (moment1 - moment0 .* 0.5.*(elevations(top_all_water_cell,1)+elevations(mixed_layer_bottom_cell+1,1))));
                    u_star_squared = 1.2 ./ density_mixed_layer.* coef_wind_drag .* wind.^2; %1.2 is density of air
                    %u_star_squared = 1.2 ./ density_mixed_layer.* 2.7e-3.*wind + 0.142e-3 .* wind.^2 + 0.0764e-3 .* wind.^3;
                    energy_wind_stir = 0.5 .* coeff_wind_stir  .* u_star_squared.^1.5 .* timestep;
                    energy_total_stir =  energy_convection + energy_wind_stir;
                    energy_avail_mix = energy_avail_mix + energy_total_stir;
                    
                    q_sqr = max(0, 2 .* energy_total_stir ./ coeff_mix_conv ./ timestep) .^(2/3);
                    density_deviation = g.*(density_water(mixed_layer_bottom_cell+1,1) - density_mixed_layer)./(0.5.*(density_mixed_layer + density_water(mixed_layer_bottom_cell+1,1)));
                    energy_uplift_cell_below_ml  = 0.5 .* (density_deviation .* sum(layerThick(top_all_water_cell:mixed_layer_bottom_cell,1)) + coeff_mix_turb .* q_sqr) .* layerThick(mixed_layer_bottom_cell+1,1);
                    
                    
                    %here the condition needs to come
                    while mixed_layer_bottom_cell < bottom_all_water_cell && energy_avail_mix >= energy_uplift_cell_below_ml
                        energy_avail_mix = energy_avail_mix - energy_uplift_cell_below_ml;
                        mixed_layer_bottom_cell = mixed_layer_bottom_cell + 1;
                        
                        %next cell
                        range = top_all_water_cell:mixed_layer_bottom_cell;
                        T_mixed_layer = sum(T(range,1) .* (layerThick(range,1)-layerThick_frozen(range,1))) ./ sum(layerThick(range,1)-layerThick_frozen(range,1));
                        T(range,1) = T_mixed_layer;
                        density_mixed_layer = water_density_UNESCO(ground, T_mixed_layer);
                        density_water(range,1) = density_mixed_layer;
                        energy_mixed_layer = sum(energy(range,1)-energy_frozen(range,1)) ./ sum(layerThick(range,1)-layerThick_frozen(range,1)) .* (layerThick(range,1)-layerThick_frozen(range,1));
                        energy(range,1) = energy_mixed_layer +  energy_frozen(range,1);
                        T(range,1) = double(energy(range,1) >0) .* T(range,1);
                        energy_frozen = double(energy<0) .* energy; %the frozen part of the energy (unfrozen part is zero)
                        layerThick_frozen = max(0, min(energy_frozen ./ (-ground.STATVAR.waterIce .* ground.CONST.L_f) .* layerThick, layerThick)); %frozen part of grid cell
                        
                        if mixed_layer_bottom_cell < bottom_all_water_cell
                            density_deviation = g.*(density_water(mixed_layer_bottom_cell+1,1) - density_mixed_layer )./(0.5.*(density_mixed_layer + density_water(mixed_layer_bottom_cell+1,1)));
                            energy_uplift_cell_below_ml = 0.5 .* (density_deviation .* sum(layerThick(range,1)-layerThick_frozen(range,1)) + coeff_mix_turb .* q_sqr) .* layerThick(mixed_layer_bottom_cell+1,1);
                        end
                        
                    end
                end
                
                %mixing from the bottom, only averaging cells
                
                %bottom layer mixing
                for mixed_layer_bottom_cell = bottom_all_water_cell:-1:top_all_water_cell+1
                    moment0 = 0;
                    moment1 = 0;
                    energy_mixing = 0;
                    density_mixed_layer = density_water(mixed_layer_bottom_cell,1);
                    
                    mixed_layer_top_cell = mixed_layer_bottom_cell;
                    
                    while mixed_layer_top_cell>1 && density_mixed_layer < density_water(mixed_layer_top_cell-1,1)
                       
                        mixed_layer_top_cell = mixed_layer_top_cell - 1;
                        range = mixed_layer_top_cell:mixed_layer_bottom_cell;
                        T_mixed_layer = sum(T(range,1) .* (layerThick(range,1)-layerThick_frozen(range,1))) ./ sum(layerThick(range,1)-layerThick_frozen(range,1));
                        density_mixed_layer = water_density_UNESCO(ground, T_mixed_layer);
                    end
                    
                    %mix the layer
                    range = mixed_layer_top_cell:mixed_layer_bottom_cell;
                    T_mixed_layer = sum(T(range,1) .* (layerThick(range,1)-layerThick_frozen(range,1))) ./ sum(layerThick(range,1)-layerThick_frozen(range,1));
                    T(range,1) = T_mixed_layer;
                    density_mixed_layer = water_density_UNESCO(ground, T_mixed_layer);
                    density_water(range,1) = density_mixed_layer;
                    energy_mixed_layer = sum(energy(range,1)-energy_frozen(range,1)) ./ sum(layerThick(range,1)-layerThick_frozen(range,1)) .* (layerThick(range,1)-layerThick_frozen(range,1));
                    energy(range,1) = energy_mixed_layer +  energy_frozen(range,1);
                    T(range,1) = double(energy(range,1) >0) .* T(range,1);
                    energy_frozen = double(energy<0) .* energy; %the frozen part of the energy (unfrozen part is zero)
                    layerThick_frozen = max(0, min(energy_frozen ./ (-ground.STATVAR.waterIce .* ground.CONST.L_f) .* layerThick, layerThick)); %frozen part of grid cell
                    
                    mixed_layer_bottom_cell = mixed_layer_top_cell;
                end
                ground.STATVAR.T = T;
                ground.STATVAR.energy = energy;
                ground.TEMP.energy_avail_mix = energy_avail_mix;
                
            end
            ground.STATVAR.density_water = density_water;
        end
        
        function ground = stratify_saline_mixed_layer(ground, tile)
            density_water = water_density_saline_UNESCO(ground, ground.STATVAR.T, ground.STATVAR.salt_conc./ground.STATVAR.area./ground.STATVAR.layerThick);
            
            top_all_water_cell = find(ground.STATVAR.energy >= 0, 1, 'first');
            bottom_all_water_cell = find(ground.STATVAR.energy >= 0, 1, 'last');
            if ~isempty(top_all_water_cell) && ~isempty(bottom_all_water_cell) && ~(top_all_water_cell==bottom_all_water_cell)
                
                T = ground.STATVAR.T;
                salt_conc = ground.STATVAR.salt_conc;
                layerThick = ground.STATVAR.layerThick;
                elevations = -cumsum([0;layerThick]);
                elevations_midpoint = -(cumsum(layerThick) - layerThick./2);
                energy = ground.STATVAR.energy;
                wind = ground.TEMP.wind .* double(energy(1,1) >= 0);
                energy_avail_mix = ground.TEMP.energy_avail_mix;
                
                coeff_mix_conv = ground.CONST.coeff_mix_conv;
                coeff_wind_stir = ground.CONST.coeff_wind_stir;
                coeff_mix_turb = ground.CONST.coeff_mix_turb;
                coef_wind_drag = ground.CONST.coef_wind_drag;
                g = ground.CONST.g;
                timestep = tile.timestep;

                mixed_layer_bottom_cell = top_all_water_cell;
                moment0 = 0;
                moment1 = 0;
                energy_mixing = 0;
                density_mixed_layer = density_water(mixed_layer_bottom_cell,1);
                
                while mixed_layer_bottom_cell < bottom_all_water_cell && density_mixed_layer >= density_water(mixed_layer_bottom_cell+1,1)
                    
                    density_x_lt = density_water(mixed_layer_bottom_cell,1) .* layerThick(mixed_layer_bottom_cell,1);
                    moment0 = moment0 + density_x_lt;
                    moment1 = moment1 + density_x_lt .* elevations_midpoint(mixed_layer_bottom_cell,1);
                    
                    mixed_layer_bottom_cell = mixed_layer_bottom_cell +1;
                    T_mixed_layer = mean(T(top_all_water_cell:mixed_layer_bottom_cell,1));
                    salt_conc_mixed_layer = sum(salt_conc(top_all_water_cell:mixed_layer_bottom_cell,1)) ./ sum(layerThick(top_all_water_cell:mixed_layer_bottom_cell,1)) ./mean(ground.STATVAR.area(top_all_water_cell:mixed_layer_bottom_cell));
                    density_mixed_layer = water_density_saline_UNESCO(ground, T_mixed_layer, salt_conc);
                    
                end
                
                %last layer
                density_x_lt = density_water(mixed_layer_bottom_cell,1) .* layerThick(mixed_layer_bottom_cell,1);
                moment0 = moment0 + density_x_lt;
                moment1 = moment1 + density_x_lt .* elevations_midpoint(mixed_layer_bottom_cell,1);
                
                %mix the layer
                T_mixed_layer = mean(T(top_all_water_cell:mixed_layer_bottom_cell,1));
                T(top_all_water_cell:mixed_layer_bottom_cell,1) = T_mixed_layer;
                salt_conc_mixed_layer = sum(salt_conc(top_all_water_cell:mixed_layer_bottom_cell,1)) ./ sum(layerThick(top_all_water_cell:mixed_layer_bottom_cell,1)) .* layerThick(top_all_water_cell:mixed_layer_bottom_cell,1);
                salt_conc(top_all_water_cell:mixed_layer_bottom_cell,1) = salt_conc_mixed_layer;
                
                density_mixed_layer = water_density_saline_UNESCO(ground, T_mixed_layer, salt_conc ./ layerThick(top_all_water_cell:mixed_layer_bottom_cell,1) ./ ground.STATVAR.area(top_all_water_cell:mixed_layer_bottom_cell,1));
                density_water(top_all_water_cell:mixed_layer_bottom_cell,1) = density_mixed_layer;
                energy_mixed_layer = sum(energy(top_all_water_cell:mixed_layer_bottom_cell,1)) ./ sum(layerThick(top_all_water_cell:mixed_layer_bottom_cell,1)) .* layerThick(top_all_water_cell:mixed_layer_bottom_cell,1);
                energy(top_all_water_cell:mixed_layer_bottom_cell,1) = energy_mixed_layer;
                
                if mixed_layer_bottom_cell < bottom_all_water_cell %not bottom reached
                    energy_convection = max(0, 0.5 .* g .* coeff_mix_conv ./density_mixed_layer .* (moment1 - moment0 .* 0.5.*(elevations(top_all_water_cell,1)+elevations(mixed_layer_bottom_cell+1,1))));
                    u_star_squared = 0.5 .* coef_wind_drag./ density_mixed_layer .* wind.^2;
                    energy_wind_stir = 0.5 .* coeff_wind_stir  .* u_star_squared.^1.5 .* timestep;
                    energy_total_stir =  energy_convection + energy_wind_stir;
                    energy_avail_mix = energy_avail_mix + energy_total_stir;
                    
                    q_sqr = max(1e-10, 2 .* energy_total_stir ./ coeff_mix_conv ./ timestep) .^(2/3);
                    density_deviation = g.*(density_water(mixed_layer_bottom_cell+1,1) - density_mixed_layer)./(0.5.*(density_mixed_layer + density_water(mixed_layer_bottom_cell+1,1)));
                    energy_uplift_cell_below_ml  = 0.5 .* (density_deviation .* sum(layerThick(top_all_water_cell:mixed_layer_bottom_cell,1)) + coeff_mix_turb .* q_sqr) .* layerThick(mixed_layer_bottom_cell+1,1);
                    
                    
                    %here the condition needs to come
                    while mixed_layer_bottom_cell < bottom_all_water_cell && energy_avail_mix >= energy_uplift_cell_below_ml
                        energy_avail_mix = energy_avail_mix - energy_uplift_cell_below_ml;
                        mixed_layer_bottom_cell = mixed_layer_bottom_cell + 1;
                        
                        %next cell
                        T_mixed_layer = mean(T(top_all_water_cell:mixed_layer_bottom_cell,1));
                        T(top_all_water_cell:mixed_layer_bottom_cell,1) = T_mixed_layer;
                        salt_conc_mixed_layer = sum(salt_conc(top_all_water_cell:mixed_layer_bottom_cell,1)) ./ sum(layerThick(top_all_water_cell:mixed_layer_bottom_cell,1)) .* layerThick(top_all_water_cell:mixed_layer_bottom_cell,1);
                        salt_conc(top_all_water_cell:mixed_layer_bottom_cell,1) = salt_conc_mixed_layer;
                        
                        energy_mixed_layer = sum(energy(top_all_water_cell:mixed_layer_bottom_cell,1)) ./ sum(layerThick(top_all_water_cell:mixed_layer_bottom_cell,1)) .* layerThick(top_all_water_cell:mixed_layer_bottom_cell,1);
                        energy(top_all_water_cell:mixed_layer_bottom_cell,1) = energy_mixed_layer;
                        density_mixed_layer = water_density_saline_UNESCO(ground, T_mixed_layer, salt_conc ./ layerThick(top_all_water_cell:mixed_layer_bottom_cell,1) ./ ground.STATVAR.area(top_all_water_cell:mixed_layer_bottom_cell,1));
                        density_water(top_all_water_cell:mixed_layer_bottom_cell,1) = density_mixed_layer;
                        
                        if mixed_layer_bottom_cell < bottom_all_water_cell
                            density_deviation = g.*(density_water(mixed_layer_bottom_cell+1,1) - density_mixed_layer )./(0.5.*(density_mixed_layer + density_water(mixed_layer_bottom_cell+1,1)));
                            energy_uplift_cell_below_ml = 0.5 .* (density_deviation .* sum(layerThick(top_all_water_cell:mixed_layer_bottom_cell,1)) + coeff_mix_turb .* q_sqr) .* layerThick(mixed_layer_bottom_cell+1,1);
                        end
                        
                    end
                end
                
                %mixing from the bottom, only averaging cells
                
                %bottom layer mixing
                for mixed_layer_bottom_cell = bottom_all_water_cell:-1:top_all_water_cell+1
                    moment0 = 0;
                    moment1 = 0;
                    energy_mixing = 0;
                    density_mixed_layer = density_water(mixed_layer_bottom_cell,1);
                    
                    mixed_layer_top_cell = mixed_layer_bottom_cell;
                    
                    while mixed_layer_top_cell>1 && density_mixed_layer < density_water(mixed_layer_top_cell-1,1)
                       
                        mixed_layer_top_cell = mixed_layer_top_cell - 1;
                        T_mixed_layer = mean(T(mixed_layer_top_cell:mixed_layer_bottom_cell,1));
                        salt_conc_mixed_layer = sum(salt_conc(mixed_layer_top_cell:mixed_layer_bottom_cell,1)) ./ sum(layerThick(mixed_layer_top_cell:mixed_layer_bottom_cell,1)) ./mean(ground.STATVAR.area(mixed_layer_top_cell:mixed_layer_bottom_cell));
                        density_mixed_layer = water_density_saline_UNESCO(ground, T_mixed_layer, salt_conc);
                    end
                    
                    %mix the layer
                    T_mixed_layer = mean(T(mixed_layer_top_cell:mixed_layer_bottom_cell,1));
                    T(mixed_layer_bottom_cell:mixed_layer_bottom_cell,1) = T_mixed_layer;
                    energy_mixed_layer = sum(energy(mixed_layer_top_cell:mixed_layer_bottom_cell,1)) ./ sum(layerThick(mixed_layer_top_cell:mixed_layer_bottom_cell,1)) .* layerThick(mixed_layer_top_cell:mixed_layer_bottom_cell,1);
                    energy(mixed_layer_top_cell:mixed_layer_bottom_cell,1) = energy_mixed_layer;
                    salt_conc_mixed_layer = sum(salt_conc(mixed_layer_top_cell:mixed_layer_bottom_cell,1)) ./ sum(layerThick(mixed_layer_top_cell:mixed_layer_bottom_cell,1)) .* layerThick(mixed_layer_top_cell:mixed_layer_bottom_cell,1);
                    salt_conc(mixed_layer_top_cell:mixed_layer_bottom_cell,1) = salt_conc_mixed_layer;
                    density_mixed_layer = water_density_saline_UNESCO(ground, T_mixed_layer, salt_conc ./ layerThick(mixed_layer_top_cell:mixed_layer_bottom_cell,1) ./ ground.STATVAR.area(mixed_layer_top_cell:mixed_layer_bottom_cell,1));
                    density_water(mixed_layer_top_cell:mixed_layer_bottom_cell,1) = density_mixed_layer;
                    
                    mixed_layer_bottom_cell = mixed_layer_top_cell;
                end
                ground.STATVAR.T = T;
                ground.STATVAR.salt_conc = salt_conc
                ground.STATVAR.energy = energy;
                ground.TEMP.energy_avail_mix = energy_avail_mix;
                
            end
            ground.STATVAR.density_water = density_water;
        end
        
        function ground = shear_mixing(ground, tile)
            wind = ground.TEMP.wind .* double(ground.STATVAR.energy(1,1) >= -0.02 .* ground.STATVAR.area(1,1) .* ground.CONST.L_f);
            ground.STATVAR.density_water = water_density_UNESCO(ground, ground.STATVAR.T);
            if wind>0                
                density_gradient = -(ground.STATVAR.density_water(1:end-2,1) - ground.STATVAR.density_water(3:end,1)) ./ (ground.STATVAR.layerThick(1:end-2,1)/2 + ground.STATVAR.layerThick(2:end-1,1) +  ground.STATVAR.layerThick(3:end,1));
                density_gradient = [ -(ground.STATVAR.density_water(1,1) - ground.STATVAR.density_water(2,1)) ./ (ground.STATVAR.layerThick(1,1)/2 + (ground.STATVAR.layerThick(2,1)/2)); density_gradient];
                density_gradient = [density_gradient; -(ground.STATVAR.density_water(end-1,1) - ground.STATVAR.density_water(end,1)) ./ (ground.STATVAR.layerThick(end-1,1)/2 + (ground.STATVAR.layerThick(end,1)/2))];
                
                N_square = ground.CONST.g ./ ground.CONST.rho_w .* density_gradient;
                
                surface_velocity = 0.028 .* wind;
                elevations = cumsum(ground.STATVAR.layerThick) - ground.STATVAR.layerThick./2;
                lake_depth = sum(ground.STATVAR.layerThick);
                velocity = surface_velocity .*( 3.* (elevations./lake_depth).^2 - 4.* elevations./lake_depth +1);
                velocity_gradient = -(velocity(1:end-2,1) - velocity(3:end,1)) ./ (ground.STATVAR.layerThick(1:end-2,1)/2 + ground.STATVAR.layerThick(2:end-1,1) +  ground.STATVAR.layerThick(3:end,1));
                velocity_gradient = [ -(velocity(1,1) - velocity(2,1)) ./ (ground.STATVAR.layerThick(1,1)/2 + (ground.STATVAR.layerThick(2,1)/2)); velocity_gradient];
                velocity_gradient = [velocity_gradient; -(velocity(end-1,1) - velocity(end,1)) ./ (ground.STATVAR.layerThick(end-1,1)/2 + (ground.STATVAR.layerThick(end,1)/2))];
                
                gradient_Richards_number = N_square ./ velocity_gradient.^2;
                
                ground.STATVAR.diffusivity_shearMixing =  1e-5 .* min(1, max(0, (1 - (gradient_Richards_number./0.7).^2).^3 ));
                
                ground.STATVAR.thermCond = ground.STATVAR.thermCond + ground.STATVAR.diffusivity_shearMixing .* ground.CONST.c_w;
%                                ground.STATVAR.thermCond = ground.STATVAR.thermCond + double(ground.STATVAR.energy >= -0.02 .* ground.STATVAR.area .* ground.CONST.L_f) .* ground.STATVAR.diffusivity_shearMixing .* ground.CONST.c_w;

            end
        end
        
        %density of water
        function density_water = water_density(ground)
            T = max(0, ground.STATVAR.T);
            %density_water = (999.83952 + 16.945176 .* T - 7.9870401e-3 .* T.^2 - 46.170461e-6 .* T.^3 + 105.56302e-9 .* T.^4 - 280.54253e-12 .* T.^5) ./ (1 + 16.897850e-3 .* T);
            %from Kell 1975, using simplified version, works well 0-30 degreeC 
            density_water = (999.83952 + 16.945176 .* T - 7.9870401e-3 .* T.^2 - 46.170461e-6 .* T.^3) ./ (1 + 16.897850e-3 .* T);
        end
        
        function density_water = water_density_UNESCO(ground, T)
%             T = max(0, ground.STATVAR.T);
            T = max(0,T);

            c1=999.842594;
            c2=6.793952e-2;
            c3=9.095290e-3;
            c4=1.001685e-4;
            c5=1.120083e-6;
            c6=6.536332e-9;
            
            t2 = T.^2;
            t3 = t2.*T;
            t4 = t3.*T;
            t5 = t4.*T;
            
            term1  =  c1;
            term2  =  c2 .* T;
            term3  = -c3 .* t2;
            term4  =  c4 .* t3;
            term5  = -c5 .* t4;
            term6  =  c6 .* t5;
            
            density_water  =  term6  + term5  + term4  + term2 + term3 + term1;
        end
        
        
        function density_water = water_density_saline_UNESCO(ground, T, salt_conc)
%             T = max(0, ground.STATVAR.T);
            T = max(0,T);
            %required salt_conc in g salt/kg water 
%             salt_conc = salt_conc .* 1e-3 .* (23+35.5)/2; %salt_conc in the model in mol/m3 = 1e-3 mol/l; assuming NaCl means half of the atoms have atomic mass 23 and half atomic mall 35.5
            salt_conc = salt_conc .* 1e-3 .* 35.5./1.1243; %salt_conc in the model in mol/m3 = 1e-3 mol/l; assuming an average ion composition of sea water, 35g/l correpsonds to 1.1243 mol/l (density of water set to 1000 for simplicity)
            
            c1=999.842594;
            c2=6.793952e-2;
            c3=9.095290e-3;
            c4=1.001685e-4;
            c5=1.120083e-6;
            c6=6.536332e-9;
            d1=8.24493e-1;
            d2=4.0899e-3;
            d3=7.6438e-5;
            d4=8.2467e-7;
            d5=5.3875e-9;
            d6=5.72466e-3;
            d7=1.0227e-4;
            d8=1.6546e-6;
            d9=4.8314e-4;
            
            t2 = T.^2;
            t3 = t2.*T;
            t4 = t3.*T;
            t5 = t4.*T;
            s2 = salt_conc.*salt_conc;
            s32 = salt_conc.^1.5;
            
            term1  =  c1;
            term2  =  c2 .* T;
            term3  = -c3 .* t2;
            term4  =  c4 .* t3;
            term5  = -c5 .* t4;
            term6  =  c6 .* t5;
            term7  =  d1;
            term8  = -d2 .* T;
            term9  =  d3 .* t2;
            term10  = -d4 .* t3;
            term11 =  d5 .* t4;
            term12 = -d6;
            term13 =  d7 .* T;
            term14 = -d8 .* t2;
            term15 =  d9;
            
            dpure  =  term6  + term5  + term4  + term2 + term3 + term1;
            csal1  = (term11 + term10  + term9  + term8 + term7) .* salt_conc;
            csal32 = (term14 + term13 + term12) .* s32;
            csal2  =  term15 .* s2;
            
            density_water = dpure + csal1 + csal32 + csal2;
        end
        
    end
end

