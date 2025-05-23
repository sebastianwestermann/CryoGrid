%========================================================================
% CryoGrid GROUND class GROUND_fcSimple_salt_seb
% heat conduction, static soil water, salt diffusion, simple freeze curve including freezing point depression due to salt
% S. Westermann, T. Ingemann-Nielsen, J. Scheer, October 2020
%========================================================================


classdef GROUND_fcSimple_salt_buoyancy_seb < SEB & HEAT_CONDUCTION & SALT & HEAT_FLUXES_LATERAL 
    
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
                
        function ground = provide_PARA(ground)
            
            ground.PARA.albedo = []; %surface albedo [-]
            ground.PARA.epsilon = []; % surface emissivity [-]
            ground.PARA.z0 = []; %roughness length [m]
            ground.PARA.rs = []; %surface resistance against evapotranspiration [sec/m] 
            ground.PARA.tortuosity = []; % tortuosity of salt diffusion [-]
            ground.PARA.tortuosity_buoyancy = [];
            ground.PARA.conductivity_function = [];

            ground.PARA.dt_max = []; %maximum possible timestep [sec]
            ground.PARA.dE_max = []; %maximum possible energy change per timestep [J/m3]
        end
        
        function ground = provide_STATVAR(ground)
            
            ground.STATVAR.upperPos = []; % upper surface elevation [m]
            ground.STATVAR.lowerPos = []; % lower surface elevation [m]
            ground.STATVAR.layerThick = []; % thickness of grid cells [m]
            ground.STATVAR.waterIce = [];  % total volume of water plus ice in a grid cell [m3]
            ground.STATVAR.mineral = []; % total volume of minerals [m3]
            ground.STATVAR.organic = []; % total volume of organics [m3]
            ground.STATVAR.energy = [];  % total internal energy[J]
            ground.STATVAR.area =[]; % area of grid cells [m2]
            
            ground.STATVAR.T = [];  % temperature [degree C]
            ground.STATVAR.water = [];  % total volume of water [m3]
            ground.STATVAR.ice = []; %total volume of ice [m3]
            ground.STATVAR.air = [];  % total volume of air [m3] - NOT USED
            ground.STATVAR.saltConc =[]; %total molar salt volume within a grid cell [mol]
            
            ground.STATVAR.thermCond = []; %thermal conductivity [W/mK]
            ground.STATVAR.deltaT = []; % freezing point depression/onset temperature of frezing for zero salt content [degree C]
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
            
            ground.CONST.k_a = [];  % thermal conductivity air
            ground.CONST.k_w = [];  % thermal conductivity water
            ground.CONST.k_i = [];  % thermal conductivity ice 
            ground.CONST.k_o = [];  % thermal conductivity organic 
            ground.CONST.k_m = [];  % thermal conductivity mineral 
            
            ground.CONST.sigma = []; % Stefan-Boltzmann constant
            ground.CONST.kappa = []; % von Karman constant
            ground.CONST.L_s = [];   %latent heat of sublimation, latent heat of evaporation handled in a dedicated function
            ground.CONST.cp = []; % specific heat capacity at constant pressure of air
            ground.CONST.g = []; % gravitational acceleration Earth surface
            ground.CONST.R = []; % universal gas constant
            ground.CONST.Tmfw = []; % freezing temperature of free water            
            ground.CONST.rho_w = []; % water density
            ground.CONST.rho_i = []; %ice density
        end
        
        function ground = convert_units(ground, tile)
                unit_converter = str2func(tile.PARA.unit_conversion_class);
                unit_converter = unit_converter();
                ground = convert_normal(unit_converter, ground, tile);
        end


        function ground = finalize_init(ground, tile)
            
            if isempty(ground.PARA.conductivity_function) || sum(isnan(ground.PARA.conductivity_function))>0
                ground.PARA.conductivity_function = 'conductivity_mixing_squares';
            end
            
            %ground = get_E_water_salt_FreezeDepress_Xice(ground); %calculate energy, water and ice contents and brine salt concentration
            ground = get_E_water_salt_freeW(ground); %calculate energy, water and ice contents and brine salt concentration

            ground = conductivity(ground); %calculate thermal conductivity
            ground = diffusivity_salt(ground); % calculate salt diffusivity 
            ground = diffusivity_salt_buoyancy(ground); % calculate salt diffusivity
            ground = water_density_saline2(ground);
            
            ground.STATVAR.Lstar = -100; %initialize with random values
            ground.STATVAR.Qh = 0;
            ground.STATVAR.Qe = 0;
           
            ground.TEMP.d_energy = ground.STATVAR.energy.*0; %set derivatives to 0
            ground.TEMP.d_salt = ground.STATVAR.energy.*0;
        end
        
        function ground = finalize_init2(ground, tile)

            ground = get_E_water_salt_freeW(ground); %calculate energy, water and ice contents and brine salt concentration
            ground = conductivity(ground); %calculate thermal conductivity
            ground = diffusivity_salt(ground); % calculate salt diffusivity 
            ground = diffusivity_salt_buoyancy(ground); % calculate salt diffusivity
            ground = water_density_saline2(ground);
        end
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)
            forcing = tile.FORCING;
            ground = surface_energy_balance(ground, forcing);
            ground = get_boundary_condition_u_ZERO_SALT(ground); %zero salt flux assumed at upper and lower boundary
        end
        
        function [ground, S_up] = penetrate_SW(ground, S_down)  %mandatory function when used in a stratigraphyw with a class that features penetration of short-wave radiation
            
            [ground, S_up] = penetrate_SW_no_transmission(ground, S_down);
        end
        
        function [ground, S_up] = penetrate_SW_PARENT(ground, S_down)  %mandatory function when used with class that features SW penetration
            [ground, S_up] = penetrate_SW_no_transmission(ground, S_down);
        end
        
        function ground = get_boundary_condition_l(ground, tile)
            forcing = tile.FORCING;
            ground.TEMP.F_lb = forcing.PARA.heatFlux_lb .* ground.STATVAR.area(end);
            ground.TEMP.d_energy(end) = ground.TEMP.d_energy(end) + ground.TEMP.F_lb;
            ground = get_boundary_condition_l_ZERO_SALT(ground);
        end
        
        function ground = get_derivatives_prognostic(ground, tile)
            ground = get_derivative_energy(ground);
            ground = get_derivative_salt_buoancy(ground);
        end
        
        function timestep = get_timestep(ground, tile) %no timestep check for salt is performed, add if necesssary
            
            timestep = get_timestep_heat_coduction(ground);
            %timestep = min(timestep, get_timestep_salt(ground));
        end
        
        function ground = advance_prognostic(ground, tile)
            timestep = tile.timestep;
            ground.STATVAR.energy = ground.STATVAR.energy + timestep .* ground.TEMP.d_energy;
            ground.STATVAR.saltConc = ground.STATVAR.saltConc + timestep .* ground.TEMP.d_salt;
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
            forcing = tile.FORCING;
            ground = L_star(ground, forcing);
        end
        
        function ground = compute_diagnostic(ground, tile)
            forcing = tile.FORCING;
            ground = get_T_water_salt_freeW(ground); % calculate temperature, water and ice contents and brine salt concentration 
            ground = conductivity(ground); %calculate thermal conductivity
            ground = diffusivity_salt(ground); % calculate salt diffusivity 
            ground = diffusivity_salt_buoyancy(ground); % calculate salt diffusivity
            ground = water_density_saline2(ground);

            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_salt = ground.STATVAR.energy.*0;
        end
        
        function ground = check_trigger(ground, tile)
            %do nothing
        end
        
        
        %-----non-mandatory functions-------
        
        function ground = surface_energy_balance(ground, forcing) %calculates the different fluxes of the surface energy balance and adds them up to get the upper boundary energy flux
            
            ground.STATVAR.Lout = (1-ground.PARA.epsilon) .* forcing.TEMP.Lin + ground.PARA.epsilon .* ground.CONST.sigma .* (ground.STATVAR.T(1)+ 273.15).^4;
            ground.STATVAR.Sout = ground.PARA.albedo .*  forcing.TEMP.Sin;
            ground.STATVAR.Qh = Q_h(ground, forcing); 
            ground.STATVAR.Qe = Q_eq(ground, forcing);
            
            ground.TEMP.F_ub = (forcing.TEMP.Sin + forcing.TEMP.Lin - ground.STATVAR.Lout - ground.STATVAR.Sout - ground.STATVAR.Qh - ground.STATVAR.Qe) .* ground.STATVAR.area(1);
            ground.TEMP.d_energy(1) = ground.TEMP.d_energy(1) + ground.TEMP.F_ub;
        end
        
        function ground = conductivity(ground)
            conductivity_function = str2func(ground.PARA.conductivity_function);
            ground = conductivity_function(ground);
        end

        
        %-----LATERAL-------------------
        
        %-------LAT3D_HEAT-------------
        function ground = lateral3D_pull_heat(ground, lateral)
            ground = lateral3D_pull_heat_simple(ground, lateral);
        end
        
        function ground = lateral3D_push_heat(ground, lateral)
            ground = lateral3D_push_heat_simple(ground, lateral);
        end
        %------------------------------
        
        %----inherited Tier 1 functions ------------
        
        function ground = get_derivative_energy(ground)
            ground = get_derivative_energy@HEAT_CONDUCTION(ground);
        end
        
        function ground = get_derivative_salt(ground)
            ground = get_derivative_salt@SALT(ground);
        end
        
        function ground = conductivity_mixing_squares(ground)
            ground = conductivity_mixing_squares@HEAT_CONDUCTION(ground);
        end
        
        function flux = Q_h(ground, forcing)
            flux = Q_h@SEB(ground, forcing);
        end
        
        function flux = Q_eq(ground, forcing)
            flux = Q_eq@SEB(ground, forcing);
        end
        
        function timestep = get_timestep_heat_coduction(ground)
            timestep = get_timestep_heat_coduction@HEAT_CONDUCTION(ground);
        end
        
        function ground = L_star(ground, forcing)
            ground = L_star@SEB(ground, forcing);
        end
        
        
        function ground = diffusivity_salt(ground)
             ground = diffusivity_salt@SALT(ground);
        end
        
        function [ground, S_up] = penetrate_SW_no_transmission(ground, S_down)
            [ground, S_up] = penetrate_SW_no_transmission@SEB(ground, S_down);
        end
        
        function gridcell_variables = get_gridcell_variables(ground)
            gridcell_variables ={'waterIce'; 'mineral'; 'organic'; 'T'; 'energy'; 'area'; 'layerThick';  'water'; 'ice'; 'saltConc'; 'salt_c_brine'; 'density_water'};
        end
        
        %-------------param file generation-----
        function ground = param_file_info(ground)
            ground = param_file_info@BASE(ground);
            %ground = provide_PARA(ground);
            
            ground.PARA.class_category = 'GROUND';
            
            %ground.PARA.options = [];
            ground.PARA.STATVAR = {'waterIce' 'mineral' 'organic' 'saltConc' 'T'};
            
            ground.PARA.default_value.albedo = {0.2};
            ground.PARA.comment.albedo = {'surface albedo [-]'};
            
            ground.PARA.default_value.epsilon = {0.99};
            ground.PARA.comment.epsilon = {'surface emissivity [-]'};
            
            ground.PARA.default_value.z0 = {0.01};
            ground.PARA.comment.z0 = {'roughness length [m]'};
            
            ground.PARA.default_value.rs = {0};
            ground.PARA.comment.rs ={'surface resistance against evapotranspiration [sec/m]'};
            
            ground.PARA.default_value.tortuosity = {1.5};
            ground.PARA.comment.tortuosity ={'tortuosity of soil for salt diffuson [-]'};
            
            ground.PARA.default_value.conductivity_function = {'conductivity_mixing_squares'};
            ground.PARA.comment.conductivity_function = {'function employed to calculate thermal conductivity, leave empty for default'};
            
            ground.PARA.default_value.dt_max = {3600};
            ground.PARA.comment.dt_max = {'maximum possible timestep [sec]'};
            
            ground.PARA.default_value.dE_max = {50000};
            ground.PARA.comment.dE_max = {'maximum possible energy change per timestep [J/m3]'};
        end
        
    end
    
end
