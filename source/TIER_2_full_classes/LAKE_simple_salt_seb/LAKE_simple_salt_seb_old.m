%========================================================================
% CryoGrid GROUND class LAKE_simple_seb
% static water body with heat conduction, free water freeze curve, surface
% energy balance
% representation of frozen water body, works in concert with
% LAKE_simple_unfrozen_seb for unfrozen water body
% S. Westermann, October 2020
%========================================================================

classdef LAKE_simple_salt_seb < SEB & HEAT_CONDUCTION & LAKE & SALT & REGRID
    
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        
        %initializes class when switching from unfrozen to frozen conditions
        function ground = initialize_from_LAKE_previous_season(ground, LAKE_simple_unfrozen)
            ground.PARA = LAKE_simple_unfrozen.PARA;
            ground.CONST = LAKE_simple_unfrozen.CONST;
            ground.STATVAR = LAKE_simple_unfrozen.STATVAR;
            %change the STATVAR
            ground = create_stratigraphy_from_STATVAR(ground);
            ground.PARA.next_season_lake_class = class(LAKE_simple_unfrozen);
        end
        
        function ground = provide_PARA(ground)
            
            ground.PARA.albedo = [];  %surface albedo [-]
            ground.PARA.epsilon = []; % surface emissivity [-]
            ground.PARA.z0 = []; %roughness length [m]
            ground.PARA.tortuosity = [];

            ground.PARA.conductivity_function = [];
            
            ground.PARA.dt_max = [];  %maximum possible timestep [sec]
            ground.PARA.dE_max = [];  %maximum possible energy change per timestep [J/m3]
            
            ground.PARA.next_season_lake_class = [];  %LAKE class that is called by check_trigger, in this case unfozen LAKE class
            % ground.PARA.next_season_lake_class = []; %class called by creation/annihilation IA class
        end
        
        function ground = provide_STATVAR(ground)
            
            ground.STATVAR.upperPos = []; % upper surface elevation [m]
            ground.STATVAR.lowerPos = []; % lower surface elevation [m]
            ground.STATVAR.layerThick = []; % thickness of grid cells [m]
            
            ground.STATVAR.waterIce = []; % total volume of water plus ice in a grid cell [m3]
            ground.STATVAR.mineral = []; % total volume of minerals [m3]
            ground.STATVAR.organic = []; % total volume of organics [m3]
            ground.STATVAR.saltConc = [];
            ground.STATVAR.energy = [];  % total internal energy [J]
            
            ground.STATVAR.T = [];  % temperature [degree C]
            ground.STATVAR.water = [];  % total volume of water [m3]
            ground.STATVAR.ice = []; %total volume of ice [m3]
            ground.STATVAR.air = [];  % total volume of air [m3] - NOT USED
            ground.STATVAR.thermCond = []; %thermal conductivity [W/mK]

            ground.STATVAR.E_f = []; %energy at which cell starts freezing
            ground.STATVAR.T_f = []; %energy at which cell starts freezing

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
            ground.CONST.R = [];
            ground.CONST.Tmfw =[];

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
            
            ground = get_E_water_salt_freeW(ground);
            ground = conductivity(ground);
            ground = diffusivity_salt_sea_ice(ground); % calculate salt diffusivity

            ground.STATVAR.Lstar = -100;
            ground.STATVAR.Qh = 0;
            ground.STATVAR.Qe = 0;
            
            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_salt = ground.STATVAR.energy.*0;
        end
        
        function ground = finalize_init2(ground, tile)

            ground = get_E_water_salt_freeW(ground);

        end
        
        
        function ground = create_stratigraphy_from_STATVAR(ground)   %create stratigraphy of ice-covered water body, with ice layer on top and cells belwo at 0 degree C 

            % layerThickSum = ground.STATVAR.layerThick;
            % ground.STATVAR.layerThick = ground.STATVAR.layerThick_store; %reassign the old grid
            % 
            % ground.STATVAR.waterIce = ground.STATVAR.waterIce ./ layerThickSum .* ground.STATVAR.layerThick;
            % ground.STATVAR.mineral = ground.STATVAR.mineral ./ layerThickSum .* ground.STATVAR.layerThick;
            % ground.STATVAR.organic = ground.STATVAR.organic ./ layerThickSum .* ground.STATVAR.layerThick;
            % ground.STATVAR.saltConc = ground.STATVAR.saltConc ./ layerThickSum .* ground.STATVAR.layerThick;
            % 
            % ground.STATVAR.air = ground.STATVAR.air ./ layerThickSum .* ground.STATVAR.layerThick;

            ground = regrid_full(ground, {'layerThick'; 'waterIce'; 'mineral'; 'organic'; 'saltConc'});

            ground.STATVAR.area = ground.STATVAR.area + ground.STATVAR.layerThick .*0;
            
            energy_below_E_f = ground.STATVAR.energy - ground.STATVAR.E_f;

            ground.STATVAR.T_f = - ground.STATVAR.saltConc ./ ground.STATVAR.waterIce .* ground.CONST.R.* ground.CONST.Tmfw.^2 ./ ground.CONST.L_f;
            ground.STATVAR.E_f =  (ground.CONST.c_w.*ground.STATVAR.waterIce + ground.CONST.c_m.*ground.STATVAR.mineral + ground.CONST.c_o.*ground.STATVAR.organic) .* ground.STATVAR.T_f; 

            ground.STATVAR.energy = ground.STATVAR.E_f;
            %assign energy to first cell
            ground.STATVAR.energy(1,1) = ground.STATVAR.energy(1,1) + energy_below_E_f;
            
            ground = get_T_water_salt_freeW(ground);
            ground = conductivity(ground);  %thermal conductivity of water
            ground = diffusivity_salt_sea_ice(ground); % calculate salt diffusivity
            
            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_salt = ground.STATVAR.energy.*0;
        end
        
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)
            forcing = tile.FORCING;
            ground = surface_energy_balance(ground, forcing);
        end
        
        function ground = get_boundary_condition_l(ground, tile)
            forcing = tile.FORCING;
            ground.TEMP.F_lb = forcing.PARA.heatFlux_lb .* ground.STATVAR.area(end);
            ground.TEMP.d_energy(end) = ground.TEMP.d_energy(end) + ground.TEMP.F_lb;
        end
        
        function [ground, S_up] = penetrate_SW(ground, S_down)  %mandatory function when used with class that features SW penetration
            [ground, S_up] = penetrate_SW_no_transmission(ground, S_down);
        end        

        function [ground, S_up] = penetrate_SW_PARENT(ground, S_down)  %mandatory function when used with class that features SW penetration
            [ground, S_up] = penetrate_SW_no_transmission(ground, S_down);
        end        
        
        function ground = get_derivatives_prognostic(ground, tile)
            if size(ground.STATVAR.energy,1)>1
                ground = get_derivative_energy(ground);
                ground = get_derivative_salt(ground);
            end
            %ground = get_derivative_salt_diffusion_convection(ground);
        end
        
        function timestep = get_timestep(ground, tile)  
            timestep = get_timestep_heat_coduction(ground);
        end
        
        function ground = advance_prognostic(ground, tile)
            timestep = tile.timestep;
            %energy
            ground.STATVAR.energy = ground.STATVAR.energy + timestep .* ground.TEMP.d_energy;
            %salt
            ground.STATVAR.saltConc = ground.STATVAR.saltConc + timestep .* ground.TEMP.d_salt;
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
            forcing = tile.FORCING;
            ground = L_star(ground, forcing);
        end
        
        function ground = compute_diagnostic(ground, tile)

            ground = get_T_water_salt_freeW(ground);
            ground.STATVAR.T_f = - ground.STATVAR.saltConc ./ ground.STATVAR.waterIce .* ground.CONST.R.* ground.CONST.Tmfw.^2 ./ ground.CONST.L_f;
            ground.STATVAR.E_f =  (ground.CONST.c_w.*ground.STATVAR.waterIce + ground.CONST.c_m.*ground.STATVAR.mineral + ground.CONST.c_o.*ground.STATVAR.organic) .* ground.STATVAR.T_f; 

            ground = move_ice_up_salt(ground);
            ground = get_T_water_salt_freeW(ground);
            ground.STATVAR.T_f = - ground.STATVAR.saltConc ./ ground.STATVAR.waterIce .* ground.CONST.R.* ground.CONST.Tmfw.^2 ./ ground.CONST.L_f;
            ground.STATVAR.E_f =  (ground.CONST.c_w.*ground.STATVAR.waterIce + ground.CONST.c_m.*ground.STATVAR.mineral + ground.CONST.c_o.*ground.STATVAR.organic) .* ground.STATVAR.T_f; 

            ground = stratify_simple_salt(ground, tile);
            ground = get_T_water_salt_freeW(ground);
            ground.STATVAR.T_f = - ground.STATVAR.saltConc ./ ground.STATVAR.waterIce .* ground.CONST.R.* ground.CONST.Tmfw.^2 ./ ground.CONST.L_f;
            ground.STATVAR.E_f =  (ground.CONST.c_w.*ground.STATVAR.waterIce + ground.CONST.c_m.*ground.STATVAR.mineral + ground.CONST.c_o.*ground.STATVAR.organic) .* ground.STATVAR.T_f; 

            ground = conductivity(ground);
            ground = diffusivity_salt_sea_ice(ground); % calculate salt diffusivity
            
            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_salt = ground.STATVAR.energy.*0;
        end
        
        %shifts to unfrozen LAKE class
        function ground = check_trigger(ground, tile)
            if sum(double(ground.STATVAR.energy<ground.STATVAR.E_f),1)==0  %all cells unfrozen
                
                ia_create_next_season_lake = get_IA_class(class(ground), ground.PARA.next_season_lake_class); %delivers IA-class that creates and initializes the next season LAKE class
                lake_next_season = create_annihilate(ia_create_next_season_lake, ground);
                
                
                %get the interaction classes from above and below
                lake_next_season.NEXT = ground.NEXT;
                lake_next_season.PREVIOUS = ground.PREVIOUS;
                lake_next_season.PREVIOUS.NEXT = lake_next_season;
                lake_next_season.NEXT.PREVIOUS = lake_next_season;
                
                %ground (CURRENT) still points to NEXT, so CURRENT.NEXT
                %will advance to the next - then it is
                %automatically handled by the garbage collection, since
                %there is no more pointer from any active variable to it
                
                %assemble new INTERACTIONS
                if ~strcmp(class(lake_next_season.PREVIOUS), 'Top')
                    ia_class = get_IA_class(class(lake_next_season.PREVIOUS), class(lake_next_season));
                    lake_next_season.IA_PREVIOUS = ia_class;
                    lake_next_season.IA_PREVIOUS.NEXT = lake_next_season;
                    lake_next_season.IA_PREVIOUS.PREVIOUS = lake_next_season.PREVIOUS;
                    lake_next_season.PREVIOUS.IA_NEXT = ia_class;
                end
                if ~strcmp(class(lake_next_season.NEXT), 'Bottom')
                    ia_class = get_IA_class(class(lake_next_season), class(lake_next_season.NEXT));
                    lake_next_season.IA_NEXT = ia_class;
                    lake_next_season.IA_NEXT.PREVIOUS = lake_next_season;
                    lake_next_season.IA_NEXT.NEXT = lake_next_season.NEXT;
                    lake_next_season.NEXT.IA_PREVIOUS = ia_class;
                end
                
            end
        end
        
        
        %-----non-mandatory functions-------
        function ground = surface_energy_balance(ground, forcing)
            ground.STATVAR.Lout = (1-ground.PARA.epsilon) .* forcing.TEMP.Lin + ground.PARA.epsilon .* ground.CONST.sigma .* (ground.STATVAR.T(1)+ 273.15).^4;
            ground.STATVAR.Sout = ground.PARA.albedo .*  forcing.TEMP.Sin;
            ground.STATVAR.Qh = Q_h(ground, forcing);
            ground.STATVAR.Qe = Q_eq_potET(ground, forcing);
            
            ground.TEMP.F_ub = (forcing.TEMP.Sin + forcing.TEMP.Lin - ground.STATVAR.Lout - ground.STATVAR.Sout - ground.STATVAR.Qh - ground.STATVAR.Qe) .* ground.STATVAR.area(1);
            ground.TEMP.d_energy(1) = ground.TEMP.d_energy(1) + ground.TEMP.F_ub;
        end
        
        function ground = conductivity(ground)
            ground = conductivity_mixing_squares(ground);
        end
        
        function gridcell_variables = get_gridcell_variables(ground)
            gridcell_variables ={'waterIce'; 'mineral'; 'organic'; 'T'; 'energy'; 'area'; 'layerThick';  'air'; 'water'; 'ice'; 'saltConc'; 'salt_c_brine'; 'E_f'; 'T_f'};
        end
        
        
        %----inherited Tier 1 functions ------------
        
        function ground = get_derivative_energy(ground)
            ground = get_derivative_energy@HEAT_CONDUCTION(ground);
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
        
        function [ground, S_up] = penetrate_SW_no_transmission(ground, S_down)
            [ground, S_up] = penetrate_SW_no_transmission@SEB(ground, S_down);
        end
        
        function timestep = get_timestep_heat_coduction(ground)
            timestep = get_timestep_heat_coduction@HEAT_CONDUCTION(ground);
        end
        
        function ground = L_star(ground, forcing)
            ground = L_star@SEB(ground, forcing);
        end
        
        function ground = get_T_water_freeW(ground)
            ground = get_T_water_freeW@HEAT_CONDUCTION(ground);
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
            
            ground.PARA.default_value.dt_max = {3600};
            ground.PARA.comment.dt_max = {'maximum possible timestep [sec]'};
            
            ground.PARA.default_value.dE_max = {50000};
            ground.PARA.comment.dE_max = {'maximum possible energy change per timestep [J/m3]'};
            
            ground.PARA.default_value.next_season_lake_class = {'LAKE_simple_unfrozen_seb'};
            ground.PARA.comment.next_season_lake_class = {'LAKE class that is called by check_trigger, in this case unfozen LAKE class'};
            
        end
    end
    
end
