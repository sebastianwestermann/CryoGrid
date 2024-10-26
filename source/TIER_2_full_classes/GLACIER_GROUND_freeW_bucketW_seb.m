%========================================================================
% CryoGrid GLACIER CLASS class GLACIER_freeW_seb
% heat conduction, free water freeze curve, surface
% energy balance
% L. Schmidt, S. Westermann, November 2021
%========================================================================

classdef GLACIER_GROUND_freeW_bucketW_seb < SEB & HEAT_CONDUCTION & SNOW & WATER_FLUXES & HEAT_FLUXES_LATERAL & WATER_FLUXES_LATERAL 

    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        

        
        function ground = provide_PARA(ground)
            
            ground.PARA.albedo_ice = [];
            ground.PARA.albedo_sediment = [];
            ground.PARA.epsilon = [];
            ground.PARA.z0 = []; %roughness length [m]
            
            ground.PARA.evaporationDepth = [];
            ground.PARA.rootDepth = [];
            ground.PARA.ratioET = [];
            
            ground.PARA.hydraulicConductivity_glacier = [];
            ground.PARA.field_capacity_glacier = [];
            ground.PARA.min_snow_age = []; %in years, snow in SNOW class is converted when older than this age 
            
            ground.PARA.target_spacing = []; %matrix
            
            ground.PARA.dt_max = [];  %maximum possible timestep [sec]
            ground.PARA.dE_max = [];  %maximum possible energy change per timestep [J/m3]

        end
        
        function ground = provide_STATVAR(ground)
            
            ground.STATVAR.upperPos = []; % upper surface elevation [m]
            ground.STATVAR.lowerPos = []; % lower surface elevation [m]
            ground.STATVAR.layerThick = []; % thickness of grid cells [m]
            ground.STATVAR.layerThick_glacier = [];
            ground.STATVAR.layerThick_sediment = [];
            
            ground.STATVAR.glacier_fraction = [];
            ground.STATVAR.ice_fraction = [];
            ground.STATVAR.waterIce = []; % total volume of water plus ice in a grid cell [m3]
            ground.STATVAR.mineral = []; % total volume of minerals [m3]
            ground.STATVAR.organic = []; % total volume of organics [m3]
            ground.STATVAR.energy = [];  % total internal energy [J]
            ground.STATVAR.field_capacity_sediment = [];
            ground.STATVAR.satHydraulicConductivity_sediment = []; %of sediment
                        
            ground.STATVAR.T = [];  % temperature [degree C]
            ground.STATVAR.water = [];  % total volume of water [m3]
            ground.STATVAR.ice = []; %total volume of ice [m3]
            ground.STATVAR.air = [];  % total volume of air [m3] - NOT USED
            ground.STATVAR.thermCond = []; %thermal conductivity [W/mK]
            ground.STATVAR.XwaterIce = [];
            ground.STATVAR.Xwater = [];
            
            ground.STATVAR.Lstar = []; %Obukhov length [m]
            ground.STATVAR.Qh = []; %sensible heat flux [W/m2]
            ground.STATVAR.Qe = []; % latent heat flux [W/m2]
            
        end
        
        function ground = provide_CONST(ground)
            
            ground.CONST.L_f = []; % volumetric latent heat of fusion, freezing
            ground.CONST.c_w = []; % volumetric heat capacity water
            ground.CONST.c_i = []; % volumetric heat capacity ice
            ground.CONST.c_o = []; % volumetric heat capacity organic
            ground.CONST.c_m = []; % volumetric heat capacity mineral
            
            ground.CONST.k_a = [];   % thermal conductivity air
            ground.CONST.k_w = [];   % thermal conductivity water
            ground.CONST.k_i = [];   % thermal conductivity ice 
            ground.CONST.k_o = [];   % thermal conductivity organic 
            ground.CONST.k_m = [];   % thermal conductivity mineral 
            
            ground.CONST.sigma = []; %Stefan-Boltzmann constant
            ground.CONST.kappa = []; % von Karman constant
            ground.CONST.L_s = [];  %latent heat of sublimation, latent heat of evaporation handled in a dedicated function
            
            ground.CONST.cp = []; %specific heat capacity at constant pressure of air
            ground.CONST.g = []; % gravitational acceleration Earth surface
            
            ground.CONST.rho_w = []; % water density
            ground.CONST.rho_i = []; %ice density
        end
        
        function ground = convert_units(ground, tile)
                
            unit_converter = str2func(tile.PARA.unit_conversion_class);
            unit_converter = unit_converter();
            ground = convert_Xice_GLACIER(unit_converter, ground, tile);
            
        end

        function ground = finalize_init(ground, tile) 

            ground.STATVAR.XwaterIce =  ground.STATVAR.waterIce .*0;
            ground.STATVAR.Xwater =  ground.STATVAR.waterIce .*0;
            
            ground.STATVAR.field_capacity = ((1-ground.STATVAR.ice_fraction) .* ground.PARA.field_capacity_glacier .* ground.STATVAR.layerThick_glacier + ...
                double(ground.STATVAR.layerThick_glacier<=0) .* ground.STATVAR.field_capacity_sediment .* ground.STATVAR.layerThick_sediment + ...
                double(ground.STATVAR.layerThick_glacier>0) .* (ground.STATVAR.layerThick_sediment - (ground.STATVAR.mineral+ ground.STATVAR.organic)./ground.STATVAR.area) ) ...
                ./ (ground.STATVAR.layerThick_sediment + ground.STATVAR.layerThick_glacier);
            
            ground = get_E_freeW(ground);
            ground.STATVAR.hydraulicConductivity =  ground.STATVAR.waterIce .*0;
            ground = hydraulic_conductivity_GLACIER_GROUND(ground); %TO DO
            
            ground.STATVAR.Lstar = -100;
            ground.STATVAR.Qh = 0;
            ground.STATVAR.Qe = 0;
            
            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_ET = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_ET_energy = ground.STATVAR.energy.*0;

            ground.TEMP.d_Xwater = ground.STATVAR.energy.*0;
            ground.TEMP.d_Xwater_energy = ground.STATVAR.energy.*0;
            
            ground.TEMP.recalculate_grid_time = -1e10; %set to a very low time
        end
        
        function ground = finalize_init2(ground, tile)

            ground = get_E_freeW(ground);
            ground = conductivity(ground);
            
            ground.TEMP.recalculate_grid_time = -1e10;
        end
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile) 
            forcing = tile.FORCING;
           
            ground = surface_energy_balance(ground, forcing);
            ground.TEMP.F_water = forcing.TEMP.rainfall ./ 1000 ./ 24 ./3600;  %possibly add water from external source here 
            ground.TEMP.T_rainWater =  max(0,forcing.TEMP.Tair);
            
            ground.TEMP.d_water(1,1) = ground.TEMP.d_water(1,1) + ground.TEMP.F_water .* ground.STATVAR.area(1,1);
            ground.TEMP.d_water_energy(1,1) = ground.TEMP.d_water_energy(1,1) + ground.TEMP.F_water .* ground.TEMP.T_rainWater .* ground.CONST.c_w .* ground.STATVAR.area(1,1);
        end
        
        function [ground, S_up] = penetrate_SW(ground, S_down)  %mandatory function when used with class that features SW penetration
            ground.TEMP.d_energy(1,1) = ground.TEMP.d_energy(1,1) + (1 - ground.STATVAR.albedo) .* sum(S_down);
            S_up = ground.STATVAR.albedo .* S_down;  % in [W/m2]
        end
        
        function [ground, S_up] = penetrate_SW_PARENT(ground, S_down)  %mandatory function when used with class that features SW penetration
            ground.TEMP.d_energy(1,1) = ground.TEMP.d_energy(1,1) + (1 - ground.STATVAR.albedo) .* sum(S_down);
            S_up = ground.STATVAR.albedo .* S_down;  % in [W/m2]
        end
        
        function ground = get_boundary_condition_l(ground, tile)
            forcing = tile.FORCING;
            ground.TEMP.F_water = 0; 
            ground.TEMP.F_lb = forcing.PARA.heatFlux_lb .* ground.STATVAR.area(end);
            ground.TEMP.d_energy(end) = ground.TEMP.d_energy(end) + ground.TEMP.F_lb;
        end
        
        function ground = get_derivatives_prognostic(ground, tile)
            ground = get_derivative_energy(ground);
            %water fluxes
            ground = get_derivative_water2(ground);
            ground = get_derivative_Xwater_GLACIER_GROUND(ground);
            %compaction of snow
            ground = compaction_firn(ground);
            
        end
        
        function timestep = get_timestep(ground, tile)
            timestep = get_timestep_heat_conduction_ice(ground);
            timestep = min(timestep, min(double(ground.TEMP.d_water <0) .* ground.STATVAR.water ./ max(1e-50, -ground.TEMP.d_water) ./10 + double(ground.TEMP.d_water >= 0) .* ground.PARA.dt_max));
%                         timestep = min(timestep, min(double(ground.TEMP.d_water <0) .* ground.STATVAR.water ./ -ground.TEMP.d_water ./10 + double(ground.TEMP.d_water >= 0) .* ground.PARA.dt_max));

        end
        
        function ground = advance_prognostic(ground, tile) 
            timestep = tile.timestep;
            %compaction
            lt_old = ground.STATVAR.layerThick_glacier;
            ice_glacier = max(0, ground.STATVAR.ice - (ground.STATVAR.layerThick_sediment .* ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic));
            ground.STATVAR.layerThick_glacier = min(ground.STATVAR.layerThick_glacier, max(ice_glacier ./ ground.STATVAR.area, ...
                ground.STATVAR.layerThick_glacier + timestep .* ground.TEMP.compact_d_D)); %reduce by compaction, but not more than zero porsity
            ground.STATVAR.ice_fraction = ice_glacier ./ ground.STATVAR.layerThick_glacier ./ ground.STATVAR.area;
            ground.STATVAR.ice_fraction(isnan(ground.STATVAR.ice_fraction) | isinf(ground.STATVAR.ice_fraction)) = 1;
            ground.STATVAR.layerThick = ground.STATVAR.layerThick + (ground.STATVAR.layerThick_glacier - lt_old);            
           
            %energy
            ground.STATVAR.energy = ground.STATVAR.energy + timestep .* (ground.TEMP.d_energy + ground.TEMP.d_water_energy +  ground.TEMP.d_Xwater_energy);
            ground.STATVAR.energy(1) = ground.STATVAR.energy(1) + timestep .* ground.TEMP.sublimation_energy;  %snowfall energy added below, when new snow layer is merged
            %mass
            ground.STATVAR.waterIce = ground.STATVAR.waterIce + timestep .* (ground.TEMP.d_water + ground.TEMP.d_Xwater);
            ground.STATVAR.layerThick = max(ground.STATVAR.layerThick_glacier + ground.STATVAR.layerThick_sediment, ground.STATVAR.layerThick + timestep .* ground.TEMP.d_Xwater ./ ground.STATVAR.area);
            
           % ground.STATVAR.layerThick(ground.TEMP.d_Xwater <0) = ground.STATVAR.layerThick(ground.TEMP.d_Xwater <0) + timestep .* ground.TEMP.d_Xwater(ground.TEMP.d_Xwater <0) ./ ground.STATVAR.area(ground.TEMP.d_Xwater <0);
            ground.STATVAR.waterIce(1) = ground.STATVAR.waterIce(1) + timestep .* ground.STATVAR.sublimation;
            ground.STATVAR.layerThick(1) = max(ground.STATVAR.layerThick_sediment(1), ground.STATVAR.layerThick(1) + timestep .* ground.STATVAR.sublimation ./ground.STATVAR.area(1,1) ./ ground.STATVAR.ice_fraction(1,1));
            ground.STATVAR.layerThick_glacier(1) = max(0, ground.STATVAR.layerThick_glacier(1) + timestep .* ground.STATVAR.sublimation ./ground.STATVAR.area(1,1) ./ ground.STATVAR.ice_fraction(1,1));
           
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
            forcing = tile.FORCING;
            ground = L_star(ground, forcing);
        end
       
        function ground = compute_diagnostic(ground, tile)     

            ground = get_T_water_freeW(ground);
            %ice in a cell melts -> decrease layerThick_glacier; water in a
            %galcier cell refreezes: increase ice_fraction
            
            ice_glacier = max(0, ground.STATVAR.ice - (ground.STATVAR.layerThick_sediment .* ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic));
            ground.STATVAR.ice_fraction = min(1, max(ground.STATVAR.ice_fraction, ice_glacier ./ ground.STATVAR.layerThick_glacier ./ ground.STATVAR.area));
            ground.STATVAR.ice_fraction(isnan(ground.STATVAR.ice_fraction) | isinf(ground.STATVAR.ice_fraction)) = 1;
            ground.STATVAR.layerThick_glacier = ice_glacier./ground.STATVAR.ice_fraction ./ ground.STATVAR.area;
            ground.STATVAR.layerThick = max(ground.STATVAR.layerThick_glacier + ground.STATVAR.layerThick_sediment, (ground.STATVAR.waterIce + ground.STATVAR.mineral + ground.STATVAR.organic)./ground.STATVAR.area);
            
            ground = modify_grid(ground, tile); %split and merge cells
            %ground = get_T_water_freeW(ground);
            
            %distribute ice over sediment, glacier and excess water
            ice_glacier = max(0, ground.STATVAR.ice - (ground.STATVAR.layerThick_sediment .* ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic));
            ground.STATVAR.ice_fraction = min(1, max(ground.STATVAR.ice_fraction, ice_glacier ./ ground.STATVAR.layerThick_glacier ./ ground.STATVAR.area));
            ground.STATVAR.ice_fraction(isnan(ground.STATVAR.ice_fraction) | isinf(ground.STATVAR.ice_fraction)) = 1;
            ground.STATVAR.layerThick_glacier = ice_glacier./ground.STATVAR.ice_fraction ./ ground.STATVAR.area;
            ground.STATVAR.layerThick = max(ground.STATVAR.layerThick_glacier + ground.STATVAR.layerThick_sediment, (ground.STATVAR.waterIce + ground.STATVAR.mineral + ground.STATVAR.organic)./ground.STATVAR.area);
            ground.STATVAR.XwaterIce = max(0, (ground.STATVAR.layerThick - ground.STATVAR.layerThick_sediment - ground.STATVAR.layerThick_glacier) .* ground.STATVAR.area .* double(ground.STATVAR.T>=0));
            ground.STATVAR.Xwater = ground.STATVAR.XwaterIce;
            
            % compute field capacity for next step
            ground.STATVAR.field_capacity = ((1-ground.STATVAR.ice_fraction) .* ground.PARA.field_capacity_glacier .* ground.STATVAR.layerThick_glacier + ...
                double(ground.STATVAR.layerThick_glacier<=0) .* ground.STATVAR.field_capacity_sediment .* ground.STATVAR.layerThick_sediment + ...
                double(ground.STATVAR.layerThick_glacier>0) .* (ground.STATVAR.layerThick_sediment - (ground.STATVAR.mineral+ ground.STATVAR.organic)./ground.STATVAR.area) ) ...
                ./ (ground.STATVAR.layerThick_sediment + ground.STATVAR.layerThick_glacier);
                       
            ground = conductivity(ground);
            ground = hydraulic_conductivity_GLACIER_GROUND(ground); %TO DO

            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_ET = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_ET_energy = ground.STATVAR.energy.*0;

            ground.TEMP.d_Xwater = ground.STATVAR.energy.*0;
            ground.TEMP.d_Xwater_energy = ground.STATVAR.energy.*0;
        end
        
        function ground = check_trigger(ground, tile)
            %do nothing
        end
        
        
        %-----non-mandatory functions-------
        function ground = surface_energy_balance(ground, forcing) %calculates the different fluxes of the surface energy balance and adds them up to get the upper boundary energy flux
            
            ground.STATVAR.albedo = (ground.STATVAR.layerThick_sediment(1,1) .* ground.PARA.albedo_sediment + ground.STATVAR.layerThick_glacier(1,1) .* ground.PARA.albedo_ice) ./ ...
                (ground.STATVAR.layerThick_sediment(1,1) + ground.STATVAR.layerThick_glacier(1,1)); 
            ground.STATVAR.albedo(isnan(ground.STATVAR.albedo) || isinf(ground.STATVAR.albedo)) = ground.PARA.albedo_ice;
            
            ground.STATVAR.Lout = (1-ground.PARA.epsilon) .* forcing.TEMP.Lin + ground.PARA.epsilon .* ground.CONST.sigma .* (ground.STATVAR.T(1)+ 273.15).^4;
            ground.STATVAR.Sout = ground.STATVAR.albedo .*  forcing.TEMP.Sin;
            ground.STATVAR.Qh = Q_h(ground, forcing);
            if ground.STATVAR.layerThick_glacier(1,1) == 0 % reacts like normal ground
                ground.STATVAR.Qe_pot = Q_eq_potET(ground, forcing);
                ground = calculateET(ground);
                ground.STATVAR.sublimation = 0;
                ground.TEMP.sublimation_energy = 0;
            else
                ground.STATVAR.Qe = Q_eq_potET_snow(ground, forcing);  %largely reacts like snow
                ground.TEMP.d_water_ET(1) = ground.TEMP.d_water_ET(1) + ground.STATVAR.evaporation;
                ground.TEMP.d_water_ET_energy(1) = ground.TEMP.d_water_ET_energy(1) + ground.TEMP.evaporation_energy;
            end
            
            ground.TEMP.F_ub = (forcing.TEMP.Sin + forcing.TEMP.Lin - ground.STATVAR.Lout - ground.STATVAR.Sout - ground.STATVAR.Qh - ground.STATVAR.Qe) .* ground.STATVAR.area(1);
            ground.TEMP.d_energy(1) = ground.TEMP.d_energy(1) + ground.TEMP.F_ub;
            
        end
        
        
        function ground = modify_grid(ground, tile)
            if tile.t > ground.TEMP.recalculate_grid_time
                ground.TEMP.recalculate_grid_time = floor(tile.t) + 10;
                ground = recalculate_target_grid(ground);
            end
            
            extensive_variables = {'layerThick'; 'layerThick_glacier'; 'layerThick_sediment'; 'waterIce'; 'XwaterIce'; 'Xwater'; 'mineral'; 'organic'; 'energy'};
            intensive_variables = {'area'; 'ice_fraction'; 'field_capacity_sediment'; 'satHydraulicConductivity_sediment'}; 
            
            split_cell = find(ground.STATVAR.layerThick_glacier + ground.STATVAR.layerThick_sediment >= 1.5 .*ground.TEMP.target_grid, 1, 'first');
            dec = 0;
            while ~isempty(split_cell)
                dec = 1;
                split_fraction = (ground.STATVAR.layerThick_glacier(split_cell,1) + ground.STATVAR.layerThick_sediment(split_cell,1) )./ ground.TEMP.target_grid(split_cell,1);
                sf1 = (split_fraction-1)./split_fraction;
                sf2 = 1./split_fraction;
                
                for j=1:size(extensive_variables,1)
                    ground.STATVAR.(extensive_variables{j,1}) = [ground.STATVAR.(extensive_variables{j,1})(1:split_cell-1,1); sf1 .*ground.STATVAR.(extensive_variables{j,1})(split_cell,1); ...
                        sf2 .* ground.STATVAR.(extensive_variables{j,1})(split_cell,1); ground.STATVAR.(extensive_variables{j,1})(split_cell+1:end,1)];
                end
                for j=1:size(intensive_variables,1)
                    ground.STATVAR.(intensive_variables{j,1}) = [ground.STATVAR.(intensive_variables{j,1})(1:split_cell-1,1); ground.STATVAR.(intensive_variables{j,1})(split_cell,1); ...
                        ground.STATVAR.(intensive_variables{j,1})(split_cell,1); ground.STATVAR.(intensive_variables{j,1})(split_cell+1:end,1)];
                end
                ground.TEMP.target_grid = [ground.TEMP.target_grid(1:split_cell-1,1); ground.TEMP.target_grid(split_cell,1); ...
                    ground.TEMP.target_grid(split_cell,1); ground.TEMP.target_grid(split_cell+1:end,1)];

                split_cell = find(ground.STATVAR.layerThick_glacier + ground.STATVAR.layerThick_sediment >= 1.5 .*ground.TEMP.target_grid, 1, 'first');
            end
            if dec
                ground = get_T_water_freeW(ground);
            end
            
            merge_cell = find(ground.STATVAR.layerThick_glacier + ground.STATVAR.layerThick_sediment < 0.5 .*ground.TEMP.target_grid, 1, 'first');
            dec = 0;
            while ~isempty(merge_cell) && size(ground.STATVAR.layerThick,1) > 1
                dec = 1;
                if merge_cell == size(ground.STATVAR.layerThick,1)
                    merge_cell = merge_cell-1;
                end
                %DO INTENSIVE CELLS FIRST, this must be done one by one, as
                %there is no common scaling variable
%                 ice_glacier = max(0, ground.STATVAR.ice(merge_cell:merge_cell+1,1) - (ground.STATVAR.layerThick_sediment(merge_cell:merge_cell+1,1) .* ground.STATVAR.area(merge_cell:merge_cell+1,1) - ground.STATVAR.mineral(merge_cell:merge_cell+1,1) - ground.STATVAR.organic(merge_cell:merge_cell+1,1)));
%                 ground.STATVAR.ice_fraction = [ground.STATVAR.ice_fraction(1:merge_cell-1,1); ...
%                     (ground.STATVAR.ice_fraction(merge_cell,1) .* ice_glacier(1,1) + ground.STATVAR.ice_fraction(merge_cell+1,1) .* ice_glacier(2,1))./ (ice_glacier(1,1)+ice_glacier(2,1)); ...
%                     ground.STATVAR.ice_fraction(merge_cell+2:end,1)];
%                 ground.STATVAR.ice_fraction(isnan(ground.STATVAR.ice_fraction) | isinf(ground.STATVAR.ice_fraction)) = 1;

                ground.STATVAR.ice_fraction = [ground.STATVAR.ice_fraction(1:merge_cell-1,1); ...
                    (ground.STATVAR.ice_fraction(merge_cell,1) .* max(1e-9, ground.STATVAR.layerThick_glacier(merge_cell,1)) + ...
                    ground.STATVAR.ice_fraction(merge_cell+1,1) .* max(1e-9, ground.STATVAR.layerThick_glacier(merge_cell+1,1)))./ (max(1e-9, ground.STATVAR.layerThick_glacier(merge_cell,1)) + max(1e-9,ground.STATVAR.layerThick_glacier(merge_cell+1,1))); ...
                    ground.STATVAR.ice_fraction(merge_cell+2:end,1)];
                ground.STATVAR.ice_fraction(isnan(ground.STATVAR.ice_fraction) | isinf(ground.STATVAR.ice_fraction)) = 1;

                ground.STATVAR.field_capacity_sediment = [ground.STATVAR.field_capacity_sediment(1:merge_cell-1,1); ...
                    (ground.STATVAR.field_capacity_sediment(merge_cell,1) .* max(1e-9, ground.STATVAR.layerThick_sediment(merge_cell,1)) + ...
                    ground.STATVAR.field_capacity_sediment(merge_cell+1,1) .* max(1e-9, ground.STATVAR.layerThick_sediment(merge_cell+1,1)))./ (max(1e-9, ground.STATVAR.layerThick_sediment(merge_cell,1)) + max(1e-9,ground.STATVAR.layerThick_sediment(merge_cell+1,1))); ...
                    ground.STATVAR.field_capacity_sediment(merge_cell+2:end,1)];
                
                %this should be redone to have proper mixing
                ground.STATVAR.satHydraulicConductivity_sediment = [ground.STATVAR.satHydraulicConductivity_sediment(1:merge_cell-1,1); ...
                    (ground.STATVAR.satHydraulicConductivity_sediment(merge_cell,1) .* max(1e-9, ground.STATVAR.layerThick_sediment(merge_cell,1)) + ...
                    ground.STATVAR.satHydraulicConductivity_sediment(merge_cell+1,1) .* max(1e-9, ground.STATVAR.layerThick_sediment(merge_cell,1)))./ (max(1e-9, ground.STATVAR.layerThick_sediment(merge_cell,1))+max(1e-9, ground.STATVAR.layerThick_sediment(merge_cell,1))); ...
                    ground.STATVAR.satHydraulicConductivity_sediment(merge_cell+2:end,1)];
                
                ground.STATVAR.area = [ground.STATVAR.area(1:merge_cell-1,1); ...
                    (ground.STATVAR.area(merge_cell,1) .* ground.STATVAR.layerThick(merge_cell,1) + ...
                    ground.STATVAR.area(merge_cell+1,1) .* ground.STATVAR.layerThick(merge_cell+1,1))./ (ground.STATVAR.layerThick(merge_cell,1)+ground.STATVAR.layerThick(merge_cell+1,1)); ...
                    ground.STATVAR.area(merge_cell+2:end,1)];
                
                ground.TEMP.target_grid = [ground.TEMP.target_grid(1:merge_cell-1,1); max(ground.TEMP.target_grid(merge_cell,1), ground.TEMP.target_grid(merge_cell+1,1)); ground.TEMP.target_grid(merge_cell+2:end,1)];
                %extensive_variables
                for j=1:size(extensive_variables,1)
                    ground.STATVAR.(extensive_variables{j,1}) = [ground.STATVAR.(extensive_variables{j,1})(1:merge_cell-1,1); ...
                        ground.STATVAR.(extensive_variables{j,1})(merge_cell,1) + ground.STATVAR.(extensive_variables{j,1})(merge_cell+1,1); ...
                        ground.STATVAR.(extensive_variables{j,1})(merge_cell+2:end,1)];
                end
                
                merge_cell = find(ground.STATVAR.layerThick_glacier + ground.STATVAR.layerThick_sediment < 0.5 .*ground.TEMP.target_grid, 1, 'first');
            end
            if dec
                ground = get_T_water_freeW(ground);
            end
        end
        
        function ground = recalculate_target_grid(ground)
            %make dependent on distance from ground surface
            depths = cumsum(ground.STATVAR.layerThick) - ground.STATVAR.layerThick ./ 2;
            ground.TEMP.target_grid = ground.STATVAR.layerThick .*0;
            for i=1:size(ground.PARA.target_spacing.upper,1)-1
                ground.TEMP.target_grid (depths >= ground.PARA.target_spacing.upper(i,1) & depths < ground.PARA.target_spacing.upper(i+1,1),1) = ground.PARA.target_spacing.spacing(i,1);
            end
            ground.TEMP.target_grid(depths >= ground.PARA.target_spacing.upper(end,1), 1) = ground.PARA.target_spacing.spacing(end,1);
            
        end
        
        
        
        function ground = conductivity(ground)
            ground = conductivity_mixing_squares(ground);
        end
        
        
        %----------LATERAL------------------
        %-----LAT_REMOVE_SURFACE_WATER-----
        function ground = lateral_push_remove_surfaceWater(ground, lateral)
            ground = lateral_push_remove_surfaceWater_GLACIER_GROUND(ground, lateral);
        end
        
        %----LAT_SEEPAGE_FACE----------
        function ground = lateral_push_remove_water_seepage(ground, lateral)
            ground = lateral_push_remove_water_seepage_GLACIER_GROUND(ground, lateral);
        end

        %----LAT_WATER_RESERVOIR------------           
        function ground = lateral_push_water_reservoir(ground, lateral)
            ground = lateral_push_water_reservoir_GLACIER_GROUND(ground, lateral);
        end
        
        %---LAT_OVERLAND_FLOW----------
        function ground = lateral_push_remove_water_overland_flow(ground, lateral)
            ground = lateral_push_water_overland_flow_GLACIER_GROUND(ground, lateral);
        end

        %---LAT3D_MASS_GLACIER------------------------------------------
        function ground = lateral3D_pull_mass_GLACIER(ground, lateral)
            
            lateral.PARENT.STATVAR.upperPos = ground.STATVAR.upperPos;
            
            lowest_ice_cell = find(ground.STATVAR.layerThick_glacier>0.01, 1, 'last');
            if isempty(lowest_ice_cell)
                lowest_ice_cell = 1;
            end
            %would be better to split off mineral part of lowest grid cell 
            lateral.PARENT.STATVAR.layerThick = ground.STATVAR.layerThick(1:lowest_ice_cell,1);
            lateral.PARENT.STATVAR.glacier_fraction = ground.STATVAR.layerThick_glacier(1:lowest_ice_cell,1) ./ ground.STATVAR.layerThick(1:lowest_ice_cell,1);
            lateral.PARENT.STATVAR.sediment_fraction = ground.STATVAR.layerThick_sediment(1:lowest_ice_cell,1) ./ ground.STATVAR.layerThick(1:lowest_ice_cell,1);
            lateral.PARENT.STATVAR.Xwater_fraction = ground.STATVAR.Xwater(1:lowest_ice_cell,1) ./ ground.STATVAR.layerThick(1:lowest_ice_cell,1)./ ground.STATVAR.area(1:lowest_ice_cell,1);            
            lateral.PARENT.STATVAR.waterIce_fraction = ground.STATVAR.waterIce(1:lowest_ice_cell,1) ./ ground.STATVAR.layerThick(1:lowest_ice_cell,1) ./ ground.STATVAR.area(1:lowest_ice_cell,1);
            lateral.PARENT.STATVAR.ice_fraction = ground.STATVAR.ice(1:lowest_ice_cell,1) ./ ground.STATVAR.layerThick(1:lowest_ice_cell,1) ./ ground.STATVAR.area(1:lowest_ice_cell,1);            
            lateral.PARENT.STATVAR.mineral_fraction = ground.STATVAR.mineral(1:lowest_ice_cell,1) ./ ground.STATVAR.layerThick(1:lowest_ice_cell,1) ./ ground.STATVAR.area(1:lowest_ice_cell,1);
            lateral.PARENT.STATVAR.organic_fraction = ground.STATVAR.organic(1:lowest_ice_cell,1) ./ ground.STATVAR.layerThick(1:lowest_ice_cell,1) ./ ground.STATVAR.area(1:lowest_ice_cell,1);
            lateral.PARENT.STATVAR.energy = ground.STATVAR.energy(1:lowest_ice_cell,1) ./ ground.STATVAR.layerThick(1:lowest_ice_cell,1) ./ ground.STATVAR.area(1:lowest_ice_cell,1);
            lateral.PARENT.STATVAR.T = ground.STATVAR.T(1:lowest_ice_cell,1);
            lateral.PARENT.STATVAR.glacier_ice_fraction = ground.STATVAR.ice_fraction(1:lowest_ice_cell,1);
            lateral.PARENT.STATVAR.field_capacity_sediment = ground.STATVAR.field_capacity_sediment(1:lowest_ice_cell,1) ; 
            lateral.PARENT.STATVAR.satHydraulicConductivity_sediment = ground.STATVAR.satHydraulicConductivity_sediment(1:lowest_ice_cell,1); 
%             lateral.PARENT.STATVAR.glacier_ice_fraction_scaled = ground.STATVAR.ice_fraction(1:lowest_ice_cell,1) .* ground.STATVAR.layerThick_glacier(1:lowest_ice_cell,1);
%             lateral.PARENT.STATVAR.field_capacity_sediment_scaled = ground.STATVAR.field_capacity_sediment(1:lowest_ice_cell,1) .* ground.STATVAR.layerThick_sediment(1:lowest_ice_cell,1); 
%             lateral.PARENT.STATVAR.satHydraulicConductivity_sediment_scaled = ground.STATVAR.satHydraulicConductivity_sediment(1:lowest_ice_cell,1).* ground.STATVAR.layerThick_sediment(1:lowest_ice_cell,1);
        end
        
        function ground = lateral3D_push_mass_GLACIER(ground, lateral)
           
            lowest_ice_cell = find(ground.STATVAR.layerThick_glacier>0.01, 1, 'last');
            if isempty(lowest_ice_cell)
                lowest_ice_cell = 1;
            end
            range = (1:lowest_ice_cell);
            layerThick_glacier = ground.STATVAR.layerThick_glacier;
            layerThick_sediment = ground.STATVAR.layerThick_sediment;
            ground.STATVAR.layerThick_glacier(range) = ground.STATVAR.layerThick_glacier(range) + lateral.PARENT.STATVAR.d_glacier ./ ground.STATVAR.area(range);
            ground.STATVAR.layerThick_sediment(range) = ground.STATVAR.layerThick_sediment(range) + lateral.PARENT.STATVAR.d_sediment ./ ground.STATVAR.area(range);
            ground.STATVAR.waterIce(range) = ground.STATVAR.waterIce(range) + lateral.PARENT.STATVAR.d_waterIce + lateral.PARENT.STATVAR.d_Xwater;
            ground.STATVAR.ice(range) = ground.STATVAR.ice(range) + lateral.PARENT.STATVAR.d_ice;
            ground.STATVAR.mineral(range) = ground.STATVAR.mineral(range) + lateral.PARENT.STATVAR.d_mineral;
            ground.STATVAR.organic(range) = ground.STATVAR.organic(range) + lateral.PARENT.STATVAR.d_organic;
            ground.STATVAR.energy(range) = ground.STATVAR.energy(range) + lateral.PARENT.STATVAR.d_energy;
            
            ice_fraction = ground.STATVAR.ice_fraction;
            ground.STATVAR.ice_fraction(range) = (ground.STATVAR.ice_fraction(range) .* layerThick_glacier(range) .* ground.STATVAR.area(range) + lateral.PARENT.STATVAR.d_glacier_ice_fraction) ./ (ground.STATVAR.layerThick_glacier(range) .* ground.STATVAR.area(range));
            ground.STATVAR.ice_fraction(isnan( ground.STATVAR.ice_fraction) | isinf(ground.STATVAR.ice_fraction)) = ice_fraction(isnan(ground.STATVAR.ice_fraction) | isinf(ground.STATVAR.ice_fraction));
            field_capacity_sediment = ground.STATVAR.field_capacity_sediment; 
            ground.STATVAR.field_capacity_sediment(range) = (ground.STATVAR.field_capacity_sediment(range) .* layerThick_sediment(range) .* ground.STATVAR.area(range) + lateral.PARENT.STATVAR.d_field_capacity_sediment) ./ (ground.STATVAR.layerThick_sediment(range).* ground.STATVAR.area(range));
            ground.STATVAR.field_capacity_sediment(isnan(ground.STATVAR.field_capacity_sediment) | isinf(ground.STATVAR.field_capacity_sediment)) = field_capacity_sediment(isnan(ground.STATVAR.field_capacity_sediment) | isinf(ground.STATVAR.field_capacity_sediment));
            satHydraulicConductivity_sediment = ground.STATVAR.satHydraulicConductivity_sediment;
            ground.STATVAR.satHydraulicConductivity_sediment(range) = (ground.STATVAR.satHydraulicConductivity_sediment(range) .* layerThick_sediment(range) .* ground.STATVAR.area(range) + lateral.PARENT.STATVAR.d_satHydraulicConductivity_sediment) ./ (ground.STATVAR.layerThick_sediment(range).* ground.STATVAR.area(range));
            ground.STATVAR.satHydraulicConductivity_sediment(isnan(ground.STATVAR.satHydraulicConductivity_sediment) | isinf(ground.STATVAR.satHydraulicConductivity_sediment)) = satHydraulicConductivity_sediment(isnan(ground.STATVAR.satHydraulicConductivity_sediment) | isinf(ground.STATVAR.satHydraulicConductivity_sediment));
        end
        
        function ground = lateral3D_pull_water_overland_flow(ground, lateral)
            ground = lateral3D_pull_water_overland_flow_XICE(ground, lateral);
            lateral.PARENT.STATVAR.depths = ground.STATVAR.upperPos; %overwritten in case it is used with a subsurface water flow class
            lateral.PARENT.STATVAR.T_water(1,1) = ground.STATVAR.T(1,1);
        end
        
        function ground = lateral3D_push_water_overland_flow(ground, lateral)
            ground.STATVAR.waterIce(1,1) = ground.STATVAR.waterIce(1,1) + lateral.PARENT.STATVAR.water_flux(1,1);
            ground.STATVAR.energy(1,1) = ground.STATVAR.energy(1,1) + lateral.PARENT.STATVAR.water_flux_energy(1,1);
        end

        
%         %-------------param file generation-----
%         function ground = param_file_info(ground)
%             ground = param_file_info@BASE(ground);
%             %ground = provide_PARA(ground);
%             
%             ground.PARA.class_category = 'GLACIER';
%             
%             %ground.PARA.options = [];
%             ground.PARA.STATVAR = {'waterIce' 'mineral' 'organic' 'T'};
%             
%             ground.PARA.default_value.albedo = {0.6};
%             ground.PARA.comment.albedo = {'surface albedo [-]'};
%             
%             ground.PARA.default_value.epsilon = {0.99};
%             ground.PARA.comment.epsilon = {'surface emissivity [-]'};
%             
%             ground.PARA.default_value.z0 = {0.01};
%             ground.PARA.comment.z0 = {'roughness length [m]'};
%             
%             ground.PARA.default_value.dt_max = {3600};
%             ground.PARA.comment.dt_max = {'maximum possible timestep [sec]'};
%             
%             ground.PARA.default_value.dE_max = {50000};
%             ground.PARA.comment.dE_max = {'maximum possible energy change per timestep [J/m3]'};
%         end

    end
    
end
