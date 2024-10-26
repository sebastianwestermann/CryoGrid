%========================================================================
% CryoGrid GROUND class LAKE_simple_bucketW_seb
% water body with heat conduction, bucket water scheme with free water evaporation, free water freeze curve, surface
% energy balance, penetration of SW radiation (single band)
% representation of frozen water body, works in concert with 
% LAKE_simple_unfrozen_bucketW_seb for unfrozen water body
% based on
% https://gmd.copernicus.org/articles/12/473/2019/gmd-12-473-2019.pdf  and
% https://hess.copernicus.org/articles/23/4969/2019/hess-23-4969-2019.pdf
% (shear mixing enhancement of thermal conductivity)
% S. Westermann, October 2020
%========================================================================

classdef LAKE_mixed_layer_bucketW_seb < SEB & HEAT_CONDUCTION & LAKE & WATER_FLUXES & REGRID & HEAT_FLUXES_LATERAL 

    
    methods

        %----mandatory functions---------------
        %----initialization--------------------
        

        
        function ground = provide_PARA(ground)
            
            ground.PARA.albedo = [];  %surface albedo [-]
            ground.PARA.epsilon = []; % surface emissivity [-]
            ground.PARA.z0 = []; %roughness length [m]
            ground.PARA.SW_extinction = []; %e-folding constant of SW extinction [1/m] 
            %spectrally resolved albedo and SW_extinction are possible,
            %must be row array - NOT TESTED!
            
            ground.PARA.dt_max = [];  %maximum possible timestep [sec]
            ground.PARA.dE_max = [];  %maximum possible energy change per timestep [J/m3]
            
            ground.PARA.threshold_water = []; %lake depth below which a trigger is called. LAKE is generally removed, depending on GROUND class below 
        end
        
        function ground = provide_STATVAR(ground)
            
            ground.STATVAR.upperPos = []; % upper surface elevation [m]
            ground.STATVAR.lowerPos = []; % lower surface elevation [m]
            ground.STATVAR.layerThick = []; % thickness of grid cells [m]
            
            ground.STATVAR.waterIce = []; % total volume of water plus ice in a grid cell [m3]
            ground.STATVAR.mineral = []; % total volume of minerals [m3]
            ground.STATVAR.organic = []; % total volume of organics [m3]
            ground.STATVAR.energy = [];  % total internal energy [J]
            
            ground.STATVAR.T = [];  % temperature [degree C]
            ground.STATVAR.water = [];  % total volume of water [m3]
            ground.STATVAR.ice = []; %total volume of ice [m3]
            ground.STATVAR.air = [];  % total volume of air [m3] - NOT USED
            ground.STATVAR.thermCond = []; %thermal conductivity [W/mK]
            
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
            
            ground.CONST.coeff_mix_conv = 0.125;
            ground.CONST.coeff_wind_stir = 0.23;
            ground.CONST.coeff_mix_turb = 0.51;
            ground.CONST.coef_wind_drag = 0.0013;
        end
        
        

        function ground = convert_units(ground, tile)
            unit_converter = str2func(tile.PARA.unit_conversion_class);
            unit_converter = unit_converter();
            ground = convert_normal(unit_converter, ground, tile);
        end
        
        function ground = finalize_init(ground, tile)
          
            ground = get_E_freeW(ground);            
            
            ground.STATVAR.Lstar = -100;
            ground.STATVAR.Qh = 0;
            ground.STATVAR.Qe = 0;
            ground.TEMP.energy_avail_mix = 0;
            
            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_energy = ground.STATVAR.energy.*0;
        end
        
        function ground = finalize_init2(ground, tile)

            ground = get_E_freeW(ground);

        end
        
        
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)
            forcing = tile.FORCING;
            ground = surface_energy_balance(ground, forcing);
            ground = get_boundary_condition_u_water_LAKE_mixed_layer(ground, forcing); %CHECK AND CHANGE HERE
            ground.TEMP.wind = forcing.TEMP.wind;
        end
        
        function [ground, S_up] = penetrate_SW(ground, S_down)  %mandatory function when used with class that features SW penetration
            [ground, S_up] = penetrate_SW_transmission_bulk(ground, S_down);
        end
        
        function [ground, S_up] = penetrate_SW_PARENT(ground, S_down)  %mandatory function when used with class that features SW penetration
            [ground, S_up] = penetrate_SW_transmission_bulk(ground, S_down);
        end
        
        function ground = get_boundary_condition_l(ground, tile)
            forcing = tile.FORCING;
            ground.TEMP.F_lb = forcing.PARA.heatFlux_lb .* ground.STATVAR.area(end);
            ground.TEMP.d_energy(end) = ground.TEMP.d_energy(end) + ground.TEMP.F_lb;
            
            ground.TEMP.F_lb_water = 0;
            ground.TEMP.F_lb_water_energy =0;
        end
        
        function ground = get_derivatives_prognostic(ground, tile)
            if size(ground.STATVAR.energy,1)>1
                ground = get_derivative_energy(ground);
            end
        end
        
        function timestep = get_timestep(ground, tile)  %could involve check for several state variables
           timestep = get_timestep_heat_coduction(ground);
        end
        
        function ground = advance_prognostic(ground, tile) 
            timestep = tile.timestep;
            %energy
            ground.STATVAR.energy = ground.STATVAR.energy + timestep .* (ground.TEMP.d_energy + ground.TEMP.d_water_energy);
            %water - no water added at top yet
            ground.STATVAR.waterIce = ground.STATVAR.waterIce + timestep .*  ground.TEMP.d_water;
            ground.STATVAR.layerThick = ground.STATVAR.layerThick + timestep .* ground.TEMP.d_water ./ ground.STATVAR.area;
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
            forcing = tile.FORCING;
            ground = L_star(ground, forcing);
        end
       
        function ground = compute_diagnostic(ground, tile)
            a=sum(ground.STATVAR.energy);
            eold = ground.STATVAR;
            ground = move_ice_up2(ground);
            if abs(a-sum(ground.STATVAR.energy))>10
               dkfellkerg 
            end
            ground = get_T_water_freeW(ground);
            ground = stratify_surface_mixed_layer(ground, tile);
            if abs(a-sum(ground.STATVAR.energy))>10
               dkfellkerg 
            end
            ground = regrid_merge(ground, {'layerThick'; 'energy'; 'waterIce'; 'mineral'; 'organic'}, {'area'; 'thermCond'}, 'waterIce');   
            if abs(a-sum(ground.STATVAR.energy))>10
               dkfellkerg 
            end
            ground = regrid_split_first_cell(ground, {'energy'; 'layerThick'; 'waterIce'; 'mineral'; 'organic'}, {'area'; 'thermCond'});
            if abs(a-sum(ground.STATVAR.energy))>10
               dkfellkerg 
            end
            
            ground = get_T_water_freeW(ground);
            ground = conductivity(ground); %molecular condution
            ground = shear_mixing(ground, tile); %enhances diffusivity in case of wind
           

            if abs(a-sum(ground.STATVAR.energy))>10
               dkfellkerg 
            end
            
            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_energy = ground.STATVAR.energy.*0;
        end
        
        
        function ground = check_trigger(ground, tile)
            
            if sum(ground.STATVAR.waterIce./ground.STATVAR.area ,1) < ground.PARA.threshold_water 
                trigger_remove_LAKE(ground.IA_NEXT, tile);
            end

        end
        
        
        %-----non-mandatory functions-------
        function ground = surface_energy_balance(ground, forcing)
            ground.STATVAR.Lout = (1-ground.PARA.epsilon) .* forcing.TEMP.Lin + ground.PARA.epsilon .* ground.CONST.sigma .* (ground.STATVAR.T(1)+ 273.15).^4;
            [ground, S_up] = penetrate_SW(ground, forcing.TEMP.Sin.* ground.STATVAR.area(1)); %distribute SW radiation
            
            ground.STATVAR.Sout = S_up ./ ground.STATVAR.area(1);
            ground.STATVAR.Qh = Q_h(ground, forcing);
            ground.STATVAR.Qe = Q_eq_potET(ground, forcing);
            %sublimation not accounted for, normally handled by the
            %snow, add if needed

            ground.TEMP.F_ub = (forcing.TEMP.Lin - ground.STATVAR.Lout - ground.STATVAR.Qh - ground.STATVAR.Qe) .* ground.STATVAR.area(1);  
            ground.TEMP.d_energy(1) = ground.TEMP.d_energy(1) + ground.TEMP.F_ub;
        end
        
        function ground = conductivity(ground)
            ground = conductivity_mixing_squares(ground);
        end
        
        
        %-----LATERAL-------------------
        %add classes for water!
        
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
            
            ground.PARA.class_category = 'LAKE';
            
            ground.PARA.STATVAR = {'waterIce' 'mineral' 'organic' 'T'};
            
            ground.PARA.default_value.albedo = {0.08};
            ground.PARA.comment.albedo = {'surface albedo [-]'};
            
            ground.PARA.default_value.epsilon = {0.99};
            ground.PARA.comment.epsilon = {'surface emissivity [-]'};
            
            ground.PARA.default_value.z0 = {0.01};
            ground.PARA.comment.z0 = {'roughness length [m]'};
            
            ground.PARA.default_value.SW_extinction = {1};
            ground.PARA.comment.SW_extinction = {'e-folding constant of SW extinction [1/m]'};
            
            ground.PARA.default_value.dt_max = {3600};
            ground.PARA.comment.dt_max = {'maximum possible timestep [sec]'};
            
            ground.PARA.default_value.dE_max = {50000};
            ground.PARA.comment.dE_max = {'maximum possible energy change per timestep [J/m3]'};
            
            ground.PARA.default_value.next_season_lake_class = {'LAKE_simple_unfrozen_bucketW_seb'};
            ground.PARA.comment.next_season_lake_class = {'LAKE class that is called by check_trigger, in this case unfozen LAKE class'};
            
            ground.PARA.default_value.threshold_water = {0.1};
            ground.PARA.comment.threshold_water = {'lake depth below which a trigger is called. LAKE is generally removed, depending on GROUND class below'};
            
        end
        
    end
    
end
