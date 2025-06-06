%========================================================================
% CryoGrid GROUND class GROUND_freezeC_RichardsEqW_seb_pressure
% heat conduction, Richards equation water scheme, freeze curve based on
% freezing=drying assumption, surface energy balance
% S. Westermann, October 2020
%========================================================================

classdef GROUND_freezeC_RichardsEqW_seb_pressure_sedimentation < SEB & HEAT_CONDUCTION & FREEZE_CURVE_KarraPainter & WATER_FLUXES & HEAT_FLUXES_LATERAL & WATER_FLUXES_LATERAL & SOIL_MECHANICS %& INITIALIZE & FREEZE_CURVE_Painter

    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
%         function ground = GROUND_freezeC_RichardsEqW_seb(index, pprovider, cprovider, forcing)  
%             ground@INITIALIZE(index, pprovider, cprovider, forcing);
%         end

        
        function ground = provide_PARA(ground)
            
            ground.PARA.albedo = [];  %surface albedo [-]
            ground.PARA.epsilon = []; % surface emissivity [-]
            ground.PARA.z0 = []; % roughness length [m] 

            ground.PARA.permeability = [];  %permeability for fluids/gases [m2]
            
            ground.PARA.dt_max = []; %maximum possible timestep [sec]
            ground.PARA.dE_max = []; %maximum possible energy change per timestep [J/m3]
            ground.PARA.dWater_max = []; %%maximum possible volumteric water content change per timestep [-] 

            ground.PARA.LUT_size_waterIce = []; %size of lookup table for the waterIce variable [-]
            ground.PARA.LUT_size_T = [];   %size of lookup table for the (temperature) T variable [-]
            ground.PARA.min_T = []; %minimum temperature for which the LUT is calculated (modeled temperatures must be above this value) [degree C]
            ground.PARA.min_waterIce = [];  %minimum waterIce value in volumetric fraction for which the LUT is calculated (modeled waterIce must be above this value) [-]
            ground.PARA.max_waterIce = []; %maximum waterIce value in volumetric fraction for which the LUT is calculated (modeled waterIce must be below this value) [-]
            ground.PARA.min_mineral_organic = [];  %maximum mineral plus organic content in volumetric fraction for which the LUT is calculated (mineral plus organic content must be below this value) [-]
            
            ground.PARA.hardBottom_cutoff = 0.03;
            ground.PARA.smoothing_factor = [];
            ground.PARA.slope_angle = []; %Slope angle [deg]
            ground.PARA.external_pressure = []; %External load from above [Pa]
            ground.PARA.sedimentation_rate = []; %Sedimentation rate [m/s]
            ground.PARA.sedimentation_id = []; %Gridcell in which sedimentation should be added
        end
        
        function ground = provide_STATVAR(ground)
            
            ground.STATVAR.upperPos = [];  % upper surface elevation [m]
            ground.STATVAR.lowerPos = [];  % lower surface elevation [m]
            ground.STATVAR.layerThick = []; % thickness of grid cells [m]
            
            ground.STATVAR.saturation = [];
            ground.STATVAR.mineral = [];   % total volume of minerals [m3]
            ground.STATVAR.organic = [];   % total volume of organics [m3]
            ground.STATVAR.Xice = []; % total volume of Xice [m3]
            ground.STATVAR.energy = [];    % total internal energy [J]
            ground.STATVAR.soil_type = []; % integer code for soil_type; 1: sand; 2: silt: 3: clay: 4: peat; 5: water (i.e. approximation of free water, very large-pore ground material).
            ground.STATVAR.initial_voidRatio = []; %=por/(1-por)  
            ground.STATVAR.compression_index = [];
            ground.STATVAR.reference_pressure = []; %100; %[Pa] correpsonds to ~10cm water column
            ground.STATVAR.permeability = [];
            
            ground.STATVAR.friction_angle_unfrozen = [];
            ground.STATVAR.friction_angle_frozen = [];
            ground.STATVAR.cohesion_unfrozen = [];
            ground.STATVAR.cohesion_frozen = [];
            
            ground.STATVAR.T = [];  % temperature [degree C]
            ground.STATVAR.water = [];  % total volume of water [m3]
            ground.STATVAR.waterPotential = []; %soil water potential [Pa]
            ground.STATVAR.ice = [];  %total volume of ice [m3]
            ground.STATVAR.air = [];  % total volume of air [m3] - NOT USED
            ground.STATVAR.thermCond = []; %thermal conductivity [W/mK]
            ground.STATVAR.hydraulicConductivity = []; % hydraulic conductivity [m/sec]
            
            ground.STATVAR.Lstar = [];  %Obukhov length [m]
            ground.STATVAR.Qh = [];     %sensible heat flux [W/m2]
            ground.STATVAR.Qe = [];     % latent heat flux [W/m2]
            
            ground.STATVAR.field_capacity = [];  %field capacity in fraction of the total volume [-]
            ground.STATVAR.excessWater = 0;  %water volume overtopping first grid cell (i.e. surface water) [m3]
        end
        
        function ground = provide_CONST(ground)
            
            ground.CONST.L_f = []; % volumetric latent heat of fusion, freezing
            ground.CONST.Tmfw = []; % freezing temperature of free water [K]
            
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
            
            ground.CONST.cp = [];  %specific heat capacity at constant pressure of air
            ground.CONST.g = [];   % gravitational acceleration Earth surface
            
            ground.CONST.R = []; %universal gas constant
            ground.CONST.molar_mass_w = []; %molar mass of water
            
            ground.CONST.rho_w = []; % water density
            ground.CONST.rho_i = []; %ice density
            ground.CONST.rho_m = []; %mineral density
            ground.CONST.rho_o = []; %organic density
            
            ground.CONST.air_pressure = 1e5; % Atmospheric pressure [Pa] -> use SI units whenever possible!
            
            %Mualem Van Genuchten model
            ground.CONST.alpha_water = [];  %alpha parameter for different soil types [m^-1]
            ground.CONST.alpha_sand = [];
            ground.CONST.alpha_silt = [];
            ground.CONST.alpha_clay = [];
            ground.CONST.alpha_peat = [];
            
            ground.CONST.n_water = [];  %n parameter for different soil types [-]
            ground.CONST.n_sand = [];
            ground.CONST.n_silt = [];
            ground.CONST.n_clay = [];
            ground.CONST.n_peat = [];
            
            ground.CONST.residual_wc_water = [];  %residual water content for different soil types [-]
            ground.CONST.residual_wc_sand = [];   %NOTE: this parameter is generally set to 0
            ground.CONST.residual_wc_silt = [];
            ground.CONST.residual_wc_clay = [];
            ground.CONST.residual_wc_peat = [];
            
            %Changed Sebastian
            %ground.CONST.density_water = 1000; % [kg/m3]
            %ground.CONST.density_ice = 1000;
            %ground.CONST.density_mineral = 2650;
            %ground.CONST.density_organic = 1300;
            %ground.CONST.air_pressure = 1e5; % Atmospheric pressure [Pa] -> use SI units whenever possible!
        end
        
        function ground = convert_units(ground, tile)
            unit_converter = str2func(tile.PARA.unit_conversion_class);
            unit_converter = unit_converter();
            %ground = convert_Xice(unit_converter, ground, tile);
            ground = convert_normal_pressure(unit_converter, ground, tile); %not with Xice, because soil is compressed without Xice
        end
         
        function ground = finalize_init(ground, tile)
            
%            ground.PARA.heatFlux_lb = tile.FORCING.PARA.heatFlux_lb;
%            ground.PARA.airT_height = tile.FORCING.PARA.airT_height;
%            ground.STATVAR.area = tile.PARA.area + ground.STATVAR.T .* 0;
            
            ground.CONST.vanGen_alpha = [ ground.CONST.alpha_sand ground.CONST.alpha_silt ground.CONST.alpha_clay ground.CONST.alpha_peat ground.CONST.alpha_water];
            ground.CONST.vanGen_n = [ ground.CONST.n_sand ground.CONST.n_silt ground.CONST.n_clay ground.CONST.n_peat ground.CONST.n_water];
            ground.CONST.vanGen_residual_wc = [ ground.CONST.residual_wc_sand ground.CONST.residual_wc_silt ground.CONST.residual_wc_clay ground.CONST.residual_wc_peat ground.CONST.residual_wc_water];
            
            %Initialization for compressed soil column
            
            %Overburden pressure
            threshold = 0.5; above_threshold = ground.STATVAR.saturation >= threshold;
            overburden_pressure_per_cell_normal = (ground.STATVAR.mineral .* ground.CONST.rho_m + ground.STATVAR.organic .* ground.CONST.rho_o + ...
                ground.STATVAR.waterIce .* ground.CONST.rho_w) ./ ground.STATVAR.area .*ground.CONST.g; %[Pa]
            overburden_pressure_per_cell_normal_boyant = (ground.STATVAR.mineral .* (ground.CONST.rho_m - ground.CONST.rho_w) + ground.STATVAR.organic .* (ground.CONST.rho_o - ground.CONST.rho_w) + ...
                + (ground.STATVAR.layerThick - ground.STATVAR.mineral - ground.STATVAR.organic - ground.STATVAR.waterIce) .* (-ground.CONST.rho_w)) ./ ground.STATVAR.area .*ground.CONST.g; %[Pa]
            overburden_pressure_per_cell = double(~above_threshold) .* overburden_pressure_per_cell_normal + ...
                double(above_threshold) .* (overburden_pressure_per_cell_normal + (ground.STATVAR.saturation-threshold)./ (1-threshold) .* (overburden_pressure_per_cell_normal_boyant - overburden_pressure_per_cell_normal));
   
            ground.STATVAR.overburden_pressure = cumsum(overburden_pressure_per_cell)-overburden_pressure_per_cell./2; %[Pa]
            ground.STATVAR.overburden_pressure = ground.STATVAR.overburden_pressure + ground.PARA.external_pressure;

            ground.STATVAR.porosity = (ground.STATVAR.initial_voidRatio - ground.STATVAR.compression_index .* log10(ground.STATVAR.overburden_pressure./ground.STATVAR.reference_pressure)) ./ ...
                (1 + ground.STATVAR.initial_voidRatio - ground.STATVAR.compression_index .* log10(ground.STATVAR.overburden_pressure./ground.STATVAR.reference_pressure));
            ground.STATVAR.porosity = min(ground.STATVAR.porosity, ground.STATVAR.initial_porosity); %do not allow higher porosities than initial void ration
            ground.STATVAR.void_ratio = ground.STATVAR.porosity ./ (1 - ground.STATVAR.porosity); %Void ratio
            ground.STATVAR.layerThick = ((ground.STATVAR.mineral + ground.STATVAR.organic) ./ (1 - ground.STATVAR.porosity)) ./ ground.STATVAR.area;
            
            %Add Xice to compressed soil column and update layerThick, porosity and void ratio
            %Attention: Xice from parameter sheet as added on top of compressed soil
            %--> this results in no excess ice if soil is compressed more than amount of Xice
            ground.STATVAR.XwaterIce = ground.STATVAR.Xice .* ground.STATVAR.initial_layerThick .* ground.STATVAR.area;
            Xice_thickness = ground.STATVAR.XwaterIce ./ ground.STATVAR.area;
            ground.STATVAR.layerThick = ground.STATVAR.layerThick + Xice_thickness;
            ground.STATVAR.Xporosity = 1 - (ground.STATVAR.mineral + ground.STATVAR.organic) ./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            ground.STATVAR.porosity = min(ground.STATVAR.Xporosity, ground.STATVAR.initial_porosity);
            ground.STATVAR.void_ratio = ground.STATVAR.porosity ./ (1 - ground.STATVAR.porosity); %Void ratio
            ground.STATVAR.Xvoid_ratio = ground.STATVAR.Xporosity ./ (1 - ground.STATVAR.Xporosity); %Void ratio
            
            %Update waterIce and saturation for compressed soil column
            waterIce = ground.STATVAR.saturation .* ground.STATVAR.Xporosity;
            ground.STATVAR.waterIce = waterIce .* ground.STATVAR.layerThick .* ground.STATVAR.area;
            ground.STATVAR.XwaterIce = max(0,ground.STATVAR.waterIce - ground.STATVAR.saturation .* ground.STATVAR.porosity .* ground.STATVAR.initial_layerThick .* ground.STATVAR.area);
            
            %Setup surfaceWaterIce
            ground.STATVAR.SurfaceWaterIce = 0;
            ground.STATVAR.SurfaceWater = 0;
            ground.STATVAR.SurfaceIce = 0;
            ground.STATVAR.SurfaceLayerThick = 0;
            ground.STATVAR.SurfaceArea = ground.STATVAR.area(1);
            ground.STATVAR.water_in_from_snow = 0;
            ground.STATVAR.water_in_from_snow_surfaceWater = 0;
            ground.STATVAR.overland_flow = 0;
            
            %Bearing capacity
            ground.STATVAR.bearing_capacity = (10.^((ground.STATVAR.initial_voidRatio - ground.STATVAR.Xporosity - ground.STATVAR.initial_voidRatio .* ground.STATVAR.Xporosity) ./ (ground.STATVAR.compression_index - ground.STATVAR.Xporosity .* ground.STATVAR.compression_index))) .* ground.STATVAR.reference_pressure;  %[Pa]
            ground.STATVAR.bearing_capacity = ground.STATVAR.bearing_capacity .* cosd(ground.PARA.slope_angle);
            ground.STATVAR.bearing_capacity = double(ground.STATVAR.T > 0) .* ground.STATVAR.bearing_capacity + double(ground.STATVAR.T <= 0) .* max(ground.STATVAR.reference_pressure, ground.STATVAR.bearing_capacity); %generate "Xice" (void ratio becomes smaller than initial void ratio)    

            %FreezeC
            ground = get_E_freezeC_Xice_pressure(ground);
            ground = conductivity(ground);
            ground = calculate_hydraulicConductivity_RichardsEq(ground);
            ground = create_LUT_freezeC(ground);
                        
            ground.STATVAR.field_capacity = (1+((ground.STATVAR.alpha./ground.CONST.g./ground.CONST.rho_w).*(ground.CONST.g.*ground.CONST.rho_w)).^ground.STATVAR.n).^(-(1-1./ground.STATVAR.n)).*ground.STATVAR.Xporosity;
            
            %Update saturation for water and ice
            ground.STATVAR.saturation = (ground.STATVAR.waterIce ./ (ground.STATVAR.layerThick .* ground.STATVAR.area)) ./ ground.STATVAR.Xporosity;

            %Calculate gravitational potential
            ground = gravitational_potential(ground);

            ground.STATVAR.Lstar = -100;
            ground.STATVAR.Qh = 0;
            ground.STATVAR.Qe = 0;
            ground.STATVAR.runoff = 0;
            %ground.STATVAR.Xwater = ground.STATVAR.waterIce .* 0;
            %ground.STATVAR.Xice = ground.STATVAR.waterIce .* 0;
            
            ground.STATVAR.year_old = []; %str2num(datestr(tile.t, 'yyyy'));
            
            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water = [0; ground.STATVAR.energy.*0]; %Extra line for surface water
            ground.TEMP.d_water_ET = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_ET_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_mineral = 0;
            ground.TEMP.d_organic = 0;
            ground.TEMP.d_air = 0;
            ground.TEMP.surface_runoff = 0;
        end
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)
            forcing = tile.FORCING;
            
            ground = surface_energy_balance(ground, forcing);
            ground = get_boundary_condition_u_RichardsEq_pressure(ground, forcing); %checked that this flux can be taken up!!
            
            id = ground.PARA.sedimentation_id;
            if sum(double(ground.STATVAR.T(1:id,1) > 0)) == id %When T < 0 in all cells above sedimentation gridcell and when there is no snow (then this function will not be used anyway)
                fraction_mineral = ground.STATVAR.mineral(id,1) ./ (ground.STATVAR.layerThick(id,1) .* ground.STATVAR.area(id,1));
                fraction_organic = ground.STATVAR.organic(id,1) ./ (ground.STATVAR.layerThick(id,1) .* ground.STATVAR.area(id,1));
                ground.TEMP.d_mineral = fraction_mineral .* ground.PARA.sedimentation_rate;
                ground.TEMP.d_organic = fraction_organic .* ground.PARA.sedimentation_rate;
                ground.TEMP.d_air = ground.PARA.sedimentation_rate - ground.TEMP.d_mineral - ground.TEMP.d_organic;
            end
        end
        
        function [ground, S_up] = penetrate_SW(ground, S_down)  %mandatory function when used with class that features SW penetration
            [ground, S_up] = penetrate_SW_no_transmission(ground, S_down);
        end

        function [ground, S_up] = penetrate_SW_PARENT(ground, S_down)  %mandatory function when used with class that features SW penetration
            [ground, S_up] = penetrate_SW_no_transmission(ground, S_down);
        end
        
        function ground = get_boundary_condition_l(ground, tile)
            forcing = tile.FORCING;
            ground.TEMP.F_lb = forcing.PARA.heatFlux_lb .* ground.STATVAR.area(end);
            ground.TEMP.d_energy(end) = ground.TEMP.d_energy(end) + ground.TEMP.F_lb;
            ground = get_boundary_condition_l_water2(ground);  %if flux not zero, check that the water flowing out is available! Not implemented here.
        end
        
        function ground = get_derivatives_prognostic(ground, tile)
            ground = get_derivative_energy(ground);
            ground = get_derivative_water_RichardsEq_pressure(ground);
        end
        
        function timestep = get_timestep(ground, tile)  
           timestep = get_timestep_heat_coduction(ground);
           timestep = min(timestep, get_timestep_water_RichardsEq_pressure(ground));
        end
        
        function ground = advance_prognostic(ground, tile) 
            timestep = tile.timestep;
            
            %Sedimentation
            id = ground.PARA.sedimentation_id;
            ground.STATVAR.mineral(id,1) = ground.STATVAR.mineral(id,1) + timestep .* ground.TEMP.d_mineral;
            ground.STATVAR.organic(id,1) = ground.STATVAR.organic(id,1) + timestep .* ground.TEMP.d_organic;
            ground.STATVAR.layerThick(id,1) = ground.STATVAR.layerThick(id,1) + (timestep .* (ground.TEMP.d_mineral + ground.TEMP.d_organic + ground.TEMP.d_air)) ./ ground.STATVAR.area(id,1);
            ground.STATVAR.Xporosity(id,1) = 1 - (ground.STATVAR.mineral(id,1) + ground.STATVAR.organic(id,1)) ./ ground.STATVAR.layerThick(id,1) ./ ground.STATVAR.area(id,1);
            ground.STATVAR.porosity(id,1) = min(ground.STATVAR.Xporosity(id,1), ground.STATVAR.initial_porosity(id,1));
            ground.STATVAR.Xvoid_ratio(id,1) = ground.STATVAR.Xporosity(id,1) ./ (1 - ground.STATVAR.Xporosity(id,1));
            ground.STATVAR.void_ratio(id,1) = ground.STATVAR.porosity(id,1) ./ (1 - ground.STATVAR.porosity(id,1));
            ground.STATVAR.saturation(id,1) = (ground.STATVAR.waterIce(id,1) ./ (ground.STATVAR.layerThick(id,1) .* ground.STATVAR.area(id,1))) ./ ground.STATVAR.Xporosity(id,1);


            prognostic_layerThick = ground.STATVAR.saturation > 1-1e-6 & (ground.TEMP.d_water(2:end,1) > 0 | (ground.TEMP.d_water(2:end,1) < 0 & ground.TEMP.no_air(2:end,1)) | (ground.TEMP.d_water(2:end,1) < 0 & ground.STATVAR.Xvoid_ratio > ground.STATVAR.initial_voidRatio));
            ground.TEMP.prognostic_layerThick = prognostic_layerThick;

            %energy
            ground.STATVAR.energy = ground.STATVAR.energy + timestep .* ground.TEMP.d_energy;
            ground.STATVAR.energy = ground.STATVAR.energy + timestep .* ground.TEMP.d_water_energy; %add energy from water advection
            
            %Update waterIce variables
            ground.STATVAR.waterIce = ground.STATVAR.waterIce + timestep .* ground.TEMP.d_water(2:end,1);
            %SurfaceWaterIce is first updated just with rainfall and evaporation to check how much water is avalable for infiltration in first ground gridcell
            %ground.STATVAR.SurfaceWaterIce = ground.STATVAR.SurfaceWaterIce + timestep .* ground.TEMP.d_rain_evap;
            %ground.STATVAR.SurfaceWater = ground.STATVAR.SurfaceWater + timestep .* ground.TEMP.d_rain_evap;
            %Check if enough surface water is available (reduce by that amount that was not available when adding it before in line above)
            if ground.STATVAR.SurfaceWater < timestep .* ground.TEMP.d_water_in_from_above(2,1)
                d_water_reduction_available_surfaceWater = timestep .* ground.TEMP.d_water_in_from_above(2,1) - ground.STATVAR.SurfaceWater;
                ground.STATVAR.waterIce(1,1) = ground.STATVAR.waterIce(1,1) - d_water_reduction_available_surfaceWater;
                ground.TEMP.d_water_in_from_above(2,1) = ground.TEMP.d_water_in_from_above(2,1) - d_water_reduction_available_surfaceWater ./ timestep; %Update variable for further calculations
            end
            %Check if enough pore space is available for infiltration from surface grid cell
            free_pore_space = ground.STATVAR.layerThick(1) .* ground.STATVAR.area(1) - ground.STATVAR.mineral(1) - ground.STATVAR.organic(1) - ground.STATVAR.waterIce(1);
            if free_pore_space < timestep .* ground.TEMP.d_water_in_from_above(2,1)
                d_water_reduction_free_pore_space = timestep .* ground.TEMP.d_water_in_from_above(2,1) - free_pore_space;
                ground.STATVAR.waterIce(1,1) = ground.STATVAR.waterIce(1,1) - d_water_reduction_free_pore_space;
                ground.TEMP.d_water_in_from_above(2,1) = ground.TEMP.d_water_in_from_above(2,1) - d_water_reduction_free_pore_space ./ timestep; %Update variable for further calculations
            end
            
            ground.STATVAR.waterIce = max(0, ground.STATVAR.waterIce);
            ground.STATVAR.XwaterIce = max(0,ground.STATVAR.waterIce - ground.STATVAR.saturation .* ground.STATVAR.porosity .* ground.STATVAR.initial_layerThick .* ground.STATVAR.area);            
            %layerThick is a prognostic variable for saturated gridcells that still carry excess pressure
            ground.STATVAR.layerThick(prognostic_layerThick) = (ground.STATVAR.waterIce(prognostic_layerThick) + ground.STATVAR.mineral(prognostic_layerThick) + ground.STATVAR.organic(prognostic_layerThick)) ./ ground.STATVAR.area(prognostic_layerThick);
            
            %Update surface water
            ground.STATVAR.SurfaceWaterIce = max(0,ground.STATVAR.SurfaceWaterIce - timestep .* ground.TEMP.d_water_in_from_above(2,1) + timestep .* ground.TEMP.d_water_in_from_below(1,1));
            ground.STATVAR.SurfaceLayerThick = ground.STATVAR.SurfaceWaterIce ./ ground.STATVAR.area(1,1);

            %ground.STATVAR.excessWater = ground.STATVAR.excessWater + timestep .* ground.TEMP.surface_runoff;
            
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
            forcing = tile.FORCING;
            ground = L_star(ground, forcing);
        end
       
        function ground = compute_diagnostic(ground, tile)
            
              %Sedimentation: make two gridcells out of one
              id = ground.PARA.sedimentation_id;
              if ground.STATVAR.layerThick(id,1) >= 2 .* ground.PARA.target_layerThick(id,1)
                  fraction_l = ground.PARA.target_layerThick(id,1) ./ ground.STATVAR.layerThick(id,1);
                  fraction_u = 1 - fraction_l;
                  
                  %Intensive variables
                  ground.STATVAR.area = [ground.STATVAR.area(1:id-1,1);ground.STATVAR.area(id,1);ground.STATVAR.area(id,1);ground.STATVAR.area(id+1:end,1)];
                  ground.STATVAR.initial_layerThick = [ground.STATVAR.initial_layerThick(1:id-1,1);ground.STATVAR.initial_layerThick(id,1);ground.STATVAR.initial_layerThick(id,1);ground.STATVAR.initial_layerThick(id+1:end,1)];
                  ground.STATVAR.initial_porosity = [ground.STATVAR.initial_porosity(1:id-1,1);ground.STATVAR.initial_porosity(id,1);ground.STATVAR.initial_porosity(id,1);ground.STATVAR.initial_porosity(id+1:end,1)];
                  ground.STATVAR.initial_voidRatio = [ground.STATVAR.initial_voidRatio(1:id-1,1);ground.STATVAR.initial_voidRatio(id,1);ground.STATVAR.initial_voidRatio(id,1);ground.STATVAR.initial_voidRatio(id+1:end,1)];
                  ground.STATVAR.porosity = [ground.STATVAR.porosity(1:id-1,1);ground.STATVAR.porosity(id,1);ground.STATVAR.porosity(id,1);ground.STATVAR.porosity(id+1:end,1)];
                  ground.STATVAR.void_ratio = [ground.STATVAR.void_ratio(1:id-1,1);ground.STATVAR.void_ratio(id,1);ground.STATVAR.void_ratio(id,1);ground.STATVAR.void_ratio(id+1:end,1)];
                  ground.STATVAR.Xporosity = [ground.STATVAR.Xporosity(1:id-1,1);ground.STATVAR.Xporosity(id,1);ground.STATVAR.Xporosity(id,1);ground.STATVAR.Xporosity(id+1:end,1)];
                  ground.STATVAR.Xvoid_ratio = [ground.STATVAR.Xvoid_ratio(1:id-1,1);ground.STATVAR.Xvoid_ratio(id,1);ground.STATVAR.Xvoid_ratio(id,1);ground.STATVAR.Xvoid_ratio(id+1:end,1)];
                  
                  if sum(isempty(ground.STATVAR.friction_angle_unfrozen)) == 0 %not available for flat terrain, only for slopes
                    ground.STATVAR.friction_angle_unfrozen = [ground.STATVAR.friction_angle_unfrozen(1:id-1,1);ground.STATVAR.friction_angle_unfrozen(id,1);ground.STATVAR.friction_angle_unfrozen(id,1);ground.STATVAR.friction_angle_unfrozen(id+1:end,1)];
                    ground.STATVAR.friction_angle_frozen = [ground.STATVAR.friction_angle_frozen(1:id-1,1);ground.STATVAR.friction_angle_frozen(id,1);ground.STATVAR.friction_angle_frozen(id,1);ground.STATVAR.friction_angle_frozen(id+1:end,1)];
                    ground.STATVAR.cohesion_unfrozen = [ground.STATVAR.cohesion_unfrozen(1:id-1,1);ground.STATVAR.cohesion_unfrozen(id,1);ground.STATVAR.cohesion_unfrozen(id,1);ground.STATVAR.cohesion_unfrozen(id+1:end,1)];
                    ground.STATVAR.cohesion_frozen = [ground.STATVAR.cohesion_frozen(1:id-1,1);ground.STATVAR.cohesion_frozen(id,1);ground.STATVAR.cohesion_frozen(id,1);ground.STATVAR.cohesion_frozen(id+1:end,1)];
                  end
                  ground.STATVAR.soil_type = [ground.STATVAR.soil_type(1:id-1,1);ground.STATVAR.soil_type(id,1);ground.STATVAR.soil_type(id,1);ground.STATVAR.soil_type(id+1:end,1)];
                  ground.STATVAR.n = [ground.STATVAR.n(1:id-1,1);ground.STATVAR.n(id,1);ground.STATVAR.n(id,1);ground.STATVAR.n(id+1:end,1)];
                  ground.STATVAR.alpha = [ground.STATVAR.alpha(1:id-1,1);ground.STATVAR.alpha(id,1);ground.STATVAR.alpha(id,1);ground.STATVAR.alpha(id+1:end,1)];
                  ground.STATVAR.compression_index = [ground.STATVAR.compression_index(1:id-1,1);ground.STATVAR.compression_index(id,1);ground.STATVAR.compression_index(id,1);ground.STATVAR.compression_index(id+1:end,1)];
                  ground.STATVAR.reference_pressure = [ground.STATVAR.reference_pressure(1:id-1,1);ground.STATVAR.reference_pressure(id,1);ground.STATVAR.reference_pressure(id,1);ground.STATVAR.reference_pressure(id+1:end,1)];
                  %ground.STATVAR.bearing_capacity = [ground.STATVAR.bearing_capacity(1:id-1,1);ground.STATVAR.bearing_capacity(id,1);ground.STATVAR.bearing_capacity(id,1);ground.STATVAR.bearing_capacity(id+1:end,1)];                  
                  ground.STATVAR.permeability = [ground.STATVAR.permeability(1:id-1,1);ground.STATVAR.permeability(id,1);ground.STATVAR.permeability(id,1);ground.STATVAR.permeability(id+1:end,1)];
                  
                  ground.STATVAR.T = [ground.STATVAR.T(1:id-1,1);ground.STATVAR.T(id,1);ground.STATVAR.T(id,1);ground.STATVAR.T(id+1:end,1)];
                  ground.STATVAR.saturation = [ground.STATVAR.saturation(1:id-1,1);ground.STATVAR.saturation(id,1);ground.STATVAR.saturation(id,1);ground.STATVAR.saturation(id+1:end,1)];
                  ground.STATVAR.hydraulicConductivity = [ground.STATVAR.hydraulicConductivity(1:id-1,1);ground.STATVAR.hydraulicConductivity(id,1);ground.STATVAR.hydraulicConductivity(id,1);ground.STATVAR.hydraulicConductivity(id+1:end,1)];
                  ground.STATVAR.field_capacity = [ground.STATVAR.field_capacity(1:id-1,1);ground.STATVAR.field_capacity(id,1);ground.STATVAR.field_capacity(id,1);ground.STATVAR.field_capacity(id+1:end,1)];
                  
                  ground.STATVAR.waterPotential = [ground.STATVAR.waterPotential(1:id-1,1);ground.STATVAR.waterPotential(id,1);ground.STATVAR.waterPotential(id,1);ground.STATVAR.waterPotential(id+1:end,1)];
                  ground.STATVAR.thermCond = [ground.STATVAR.thermCond(1:id-1,1);ground.STATVAR.thermCond(id,1);ground.STATVAR.thermCond(id,1);ground.STATVAR.thermCond(id+1:end,1)];
                  ground.STATVAR.viscosity_water = [ground.STATVAR.viscosity_water(1:id-1,1);ground.STATVAR.viscosity_water(id,1);ground.STATVAR.viscosity_water(id,1);ground.STATVAR.viscosity_water(id+1:end,1)];

                  %Extensive variable
                  layerThick_l = ground.PARA.target_layerThick(id,1);
                  layerThick_u = ground.STATVAR.layerThick(id,1) - layerThick_l;
                  ground.STATVAR.layerThick = [ground.STATVAR.layerThick(1:id-1,1);layerThick_u;layerThick_l;ground.STATVAR.layerThick(id+1:end,1)];
                  
                  mineral_l = fraction_l .* ground.STATVAR.mineral(id,1);
                  mineral_u = ground.STATVAR.mineral(id,1) - mineral_l;
                  ground.STATVAR.mineral = [ground.STATVAR.mineral(1:id-1,1);mineral_u;mineral_l;ground.STATVAR.mineral(id+1:end,1)];
                  
                  organic_l = fraction_l .* ground.STATVAR.organic(id,1);
                  organic_u = ground.STATVAR.organic(id,1) - organic_l;
                  ground.STATVAR.organic = [ground.STATVAR.organic(1:id-1,1);organic_u;organic_l;ground.STATVAR.organic(id+1:end,1)];
                  
                  waterIce_l = fraction_l .* ground.STATVAR.waterIce(id,1);
                  waterIce_u = ground.STATVAR.waterIce(id,1) - waterIce_l;
                  ground.STATVAR.waterIce = [ground.STATVAR.waterIce(1:id-1,1);waterIce_u;waterIce_l;ground.STATVAR.waterIce(id+1:end,1)];
                  
                  XwaterIce_l = fraction_l .* ground.STATVAR.XwaterIce(id,1);
                  XwaterIce_u = ground.STATVAR.XwaterIce(id,1) - XwaterIce_l;
                  ground.STATVAR.XwaterIce = [ground.STATVAR.XwaterIce(1:id-1,1);XwaterIce_u;XwaterIce_l;ground.STATVAR.XwaterIce(id+1:end,1)];
                                  
                  water_l = fraction_l .* ground.STATVAR.water(id,1);
                  water_u = ground.STATVAR.water(id,1) - water_l;
                  ground.STATVAR.water = [ground.STATVAR.water(1:id-1,1);water_u;water_l;ground.STATVAR.water(id+1:end,1)];
                  
                  Xwater_l = fraction_l .* ground.STATVAR.Xwater(id,1);
                  Xwater_u = ground.STATVAR.Xwater(id,1) - Xwater_l;
                  ground.STATVAR.Xwater = [ground.STATVAR.Xwater(1:id-1,1);Xwater_u;Xwater_l;ground.STATVAR.Xwater(id+1:end,1)];
                  
                  ice_l = fraction_l .* ground.STATVAR.ice(id,1);
                  ice_u = ground.STATVAR.ice(id,1) - ice_l;
                  ground.STATVAR.ice = [ground.STATVAR.ice(1:id-1,1);ice_u;ice_l;ground.STATVAR.ice(id+1:end,1)];
                  
                  Xice_l = fraction_l .* ground.STATVAR.Xice(id,1);
                  Xice_u = ground.STATVAR.Xice(id,1) - Xice_l;
                  ground.STATVAR.Xice = [ground.STATVAR.Xice(1:id-1,1);Xice_u;Xice_l;ground.STATVAR.Xice(id+1:end,1)];
                  
                  air_l = fraction_l .* ground.STATVAR.air(id,1);
                  air_u = ground.STATVAR.air(id,1) - air_l;
                  ground.STATVAR.air = [ground.STATVAR.air(1:id-1,1);air_u;air_l;ground.STATVAR.air(id+1:end,1)];
                  
                  energy_l = fraction_l .* ground.STATVAR.energy(id,1);
                  energy_u = ground.STATVAR.energy(id,1) - energy_l;
                  ground.STATVAR.energy = [ground.STATVAR.energy(1:id-1,1);energy_u;energy_l;ground.STATVAR.energy(id+1:end,1)];
                    
                  %bearing_capacity
                  %overburden_pressure
                  %gravitational_potential          
              end
            
            ground = soil_mechanics(ground, tile);
            if ground.PARA.slope_angle > 0
                ground = calculate_FOS(ground, tile);
            end
            
            ground = get_T_water_freezeC_Xice_pressure(ground);

            ground = gravitational_potential(ground);
            ground = conductivity(ground);
            ground = calculate_hydraulicConductivity_RichardsEq(ground);

            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water = [0; ground.STATVAR.energy.*0]; %Extra line for surface water
            ground.TEMP.d_water_ET = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_ET_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_mineral = 0;
            ground.TEMP.d_organic = 0;
            ground.TEMP.d_air = 0;
            ground.TEMP.surface_runoff = 0;
            
        end
        
         function ground = check_trigger(ground, tile)
%             % Every year, the heat reservoir is 10% further away
%             if isempty(ground.STATVAR.year_old) || str2num(datestr(tile.t, 'yyyy'))>ground.STATVAR.year_old
%                 ground.STATVAR.year_old = str2num(datestr(tile.t, 'yyyy'));
%                 tile.LATERAL.IA_CLASSES{2,1}.PARA.distance_heatReservoir = tile.LATERAL.IA_CLASSES{2,1}.PARA.distance_heatReservoir .*1.1;
%             end 
         end
        
        
        %-----non-mandatory functions-------
        
        function ground = surface_energy_balance(ground, forcing)
            ground.STATVAR.Lout = (1-ground.PARA.epsilon) .* forcing.TEMP.Lin + ground.PARA.epsilon .* ground.CONST.sigma .* (ground.STATVAR.T(1)+ 273.15).^4;
            ground.STATVAR.Sout = ground.PARA.albedo .*  forcing.TEMP.Sin;
            ground.STATVAR.Qh = Q_h(ground, forcing);
            ground = Q_evap_CLM4_5(ground, forcing);
            
            ground.TEMP.F_ub = (forcing.TEMP.Sin + forcing.TEMP.Lin - ground.STATVAR.Lout - ground.STATVAR.Sout - ground.STATVAR.Qh - ground.STATVAR.Qe) .* ground.STATVAR.area(1);
            ground.TEMP.d_energy(1) = ground.TEMP.d_energy(1) + ground.TEMP.F_ub;

            %             %water -> evaporation
            %             ground.TEMP.d_water_ET(1,1) = ground.TEMP.d_water_ET(1,1) -  ground.STATVAR.evap.* ground.STATVAR.area(1,1); %in m3 water per sec, put everything in uppermost grid cell
            %             ground.TEMP.d_water_ET_energy(1,1) = ground.TEMP.d_water_ET_energy(1,1) -  ground.STATVAR.evap_energy.* ground.STATVAR.area(1,1);
            %water -> evaporation
            ground.TEMP.d_water_ET(1,1) = ground.TEMP.d_water_ET(1,1) +  ground.STATVAR.evaporation.* ground.STATVAR.area(1,1); %in m3 water per sec, put everything in uppermost grid cell
            ground.TEMP.d_water_ET_energy(1,1) = ground.TEMP.d_water_ET_energy(1,1) +  ground.TEMP.evaporation_energy.* ground.STATVAR.area(1,1);

        end
        
        function ground = conductivity(ground)
            ground = thermalConductivity_CLM4_5(ground);
            %ground = conductivity_mixing_squares(ground);
        end
        
        %-----LATERAL-------------------
        
        %-----LAT_REMOVE_SURFACE_WATER-----
        function ground = lateral_push_remove_surfaceWater(ground, lateral)
            ground = lateral_push_remove_surfaceWater_simple(ground, lateral);
        end
        
        %-----LAT_REMOVE_SUBSURFACE_WATER-----        
        function ground = lateral_push_remove_subsurfaceWater(ground, lateral)
            ground = lateral_push_remove_subsurfaceWater_simple(ground, lateral);
        end
        
        %----LAT_SEEPAGE_FACE----------            
        function ground = lateral_push_remove_water_seepage(ground, lateral)
            ground = lateral_push_remove_water_seepage_simple(ground, lateral);
        end
        
        %----LAT_WATER_RESERVOIR------------          
        function ground = lateral_push_water_reservoir(ground, lateral)
            ground = lateral_push_water_reservoir_RichardsEq_pressure(ground, lateral);
        end
        
        %---LAT_OVERLAND_FLOW----------
        function ground = lateral_push_remove_water_overland_flow(ground, lateral)
            ground = lateral_push_water_overland_flow_RichardsEq_pressure(ground, lateral);
        end
        
        %----LAT_HEAT------------
        function ground = lateral_push_heat(ground, lateral)
            ground = lateral_push_heat_simple(ground, lateral);
        end

        %----LAT3D_WATER_UNCONFINED_AQUIFER------------         
        function ground = lateral3D_pull_water_unconfined_aquifer(ground, lateral)
            ground = lateral3D_pull_water_unconfined_aquifer_simple(ground, lateral);
        end
        
        function ground = lateral3D_push_water_unconfined_aquifer(ground, lateral)
            ground = lateral3D_push_water_unconfined_aquifer_simple(ground, lateral);
        end
        
        function [saturated_next, hardBottom_next] = get_saturated_hardBottom_first_cell(ground, lateral)
            [saturated_next, hardBottom_next] = get_saturated_hardBottom_first_cell_simple(ground, lateral);
        end
        
        %LAT3D_WATER_RESERVOIR and LAT3D_WATER_SEEPAGE_FACE do not require specific functions
         
        %-------LAT3D_HEAT-------------        
        function ground = lateral3D_pull_heat(ground, lateral)
            ground = lateral3D_pull_heat_simple(ground, lateral);
        end
        
        function ground = lateral3D_push_heat(ground, lateral)
            ground = lateral3D_push_heat_simple(ground, lateral);
        end
        
                        %-------------param file generation-----
        function ground = param_file_info(ground)
            ground = param_file_info@BASE(ground);
            
            ground.PARA.class_category = 'GROUND';
            
            ground.PARA.STATVAR = {'waterIce' 'mineral' 'organic' 'soil_type' 'field_capacity' 'permeability' 'T'}; %add further variables here
            
            ground.PARA.default_value.albedo = {0.2};
            ground.PARA.comment.albedo = {'surface albedo [-]'};
            
            ground.PARA.default_value.epsilon = {0.99};
            ground.PARA.comment.epsilon = {'surface emissivity [-]'};
            
            ground.PARA.default_value.z0 = {0.01};
            ground.PARA.comment.z0 = {'roughness length [m]'};
            
            ground.PARA.default_value.rootDepth = {0.1};
            ground.PARA.comment.rootDepth = {'e-folding depth of transpiration reduction with depth [m]'};
            
            ground.PARA.default_value.evaporationDepth = {0.1};
            ground.PARA.comment.evaporationDepth = {'e-folding constant of evaporation reduction reduction with depth [m]'};
            
            ground.PARA.default_value.ratioET = {0.5};
            ground.PARA.comment.ratioET = {'fraction of transpiration of total evapotranspiration [-]'};
            
            ground.PARA.default_value.conductivity_function = {''};
            ground.PARA.comment.conductivity_function = {'function employed to calculate thermal conductivity, leave empty for default'};
            
            ground.PARA.default_value.dt_max = {3600};
            ground.PARA.comment.dt_max = {'maximum possible timestep [sec]'};
            
            ground.PARA.default_value.dE_max = {50000};
            ground.PARA.comment.dE_max = {'maximum possible energy change per timestep [J/m3]'};
            
            ground.PARA.default_value.dWater_max = {0.005};
            ground.PARA.comment.dWater_max = {'maximum possible volumteric water content change per timestep [-]'};
            
            ground.PARA.default_value.LUT_size_waterIce = {1000};
            ground.PARA.comment.LUT_size_waterIce = {'size of lookup table for the waterIce variable [-]'};
            
            ground.PARA.default_value.LUT_size_T = {1000};
            ground.PARA.comment.LUT_size_T = {'size of lookup table for the (temperature) T variable [-]'};
            
            ground.PARA.default_value.min_T = {-50};
            ground.PARA.comment.min_T = {'minimum temperature for which the LUT is calculated (modeled temperatures must be above this value) [degree C]'};
            
            ground.PARA.default_value.min_waterIce = {0.05};
            ground.PARA.comment.min_waterIce = {'minimum waterIce value in volumetric fraction for which the LUT is calculated (modeled waterIce must be above this value) [-]'};
            
            ground.PARA.default_value.max_waterIce = {0.97};
            ground.PARA.comment.max_waterIce = {'maximum waterIce value in volumetric fraction for which the LUT is calculated (modeled waterIce must be below this value) [-]'};
            
            ground.PARA.default_value.min_mineral_organic = {0.03};
            ground.PARA.comment.min_mineral_organic = {'maximum mineral plus organic content in volumetric fraction for which the LUT is calculated (mineral plus organic content must be below this value) [-]'};
        
            ground.PARA.default_value.hardBottom_cutoff = {0.03};
            ground.PARA.comment.hardBottom_cutoff = {'ground considered impermeable for water when water content is less than cutoff'};
            
            ground.PARA.default_value.smoothing_factor = {100};
            ground.PARA.comment.smoothing_factor = {'adjustment timescale of bearing capacity'};

        end
        
    end
    
end
