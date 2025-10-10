%========================================================================
% CryoGrid GROUND class GROUND_freeW_heatPump_seb
% heat conduction, free water freeze curve, constant water + ice water balance, 
% surface energy balance
% S. Westermann, Aug 2025
%========================================================================

classdef GROUND_freeW_heatPump_seb < SEB & HEAT_CONDUCTION & HEATPUMP & HEAT_FLUXES_LATERAL & ADJUST_STRATIGRAPHY
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        
        function ground = provide_PARA(ground)
            
            ground.PARA.albedo = []; %surface albedo [-]
            ground.PARA.epsilon = []; % surface emissivity [-]
            ground.PARA.z0 = []; %roughness length [m]
            ground.PARA.rs = []; %surface resistance against evapotranspiration [secm] 

            ground.PARA.number_of_pipes = []; %total number, no matter the geometry
            ground.PARA.number_of_circular_cells = [];
            ground.PARA.distance_between_pipes = [];
            ground.PARA.pipe_diameter = [];
            ground.PARA.heatCapacity_fluid = 4.2e6;
            ground.PARA.fluid_speed = []; %>0 flow top down, <0 flow bottom up; %flow speed is flow speed in each indiuvidual pipe
            ground.PARA.incoming_fluid_T = [];
            ground.PARA.max_power = 5000; %[W]

            ground.PARA.conductivity_function = [];
            ground.PARA.dt_max = []; %maximum possible timestep [sec]
            ground.PARA.dE_max = []; %maximum possible energy change per timestep [J/m3]
        end
        
        function ground = provide_STATVAR(ground)
            
            ground.STATVAR.upperPos = []; % upper surface elevation [m]
            ground.STATVAR.lowerPos = []; % lower surface elevation [m]
            ground.STATVAR.layerThick = [];  % thickness of grid cells [m]
            
            ground.STATVAR.waterIce = []; % total volume of water plus ice in a grid cell [m3]
            ground.STATVAR.mineral = [];  % total volume of minerals [m3]
            ground.STATVAR.organic = [];  %total volume of ice [m3]
            ground.STATVAR.cooling_fluid = [];
            ground.STATVAR.heatPump_active = [];
            ground.STATVAR.energy = [];  %total molar salt volume within a grid cell [mol]
            
            ground.STATVAR.T = [];   % temperature [degree C]
            ground.STATVAR.water = [];  % total volume of water [m3]
            ground.STATVAR.ice = [];    %total volume of ice [m3]
            ground.STATVAR.thermCond = []; %thermal conductivity [W/mK]
            ground.STATVAR.fluid_speed = [];
            ground.STATVAR.incoming_fluid_T = [];
            ground.STATVAR.outgoing_fluid_T = [];
            ground.STATVAR.area = [];
            ground.STATVAR.contact_length = [];
            
            ground.STATVAR.Lstar = [];  %Obukhov length [m]
            ground.STATVAR.Qh = []; %sensible heat flux [W/m2]
            ground.STATVAR.Qe = [];  % latent heat flux [W/m2]
        end
    
        function ground = provide_CONST(ground)
            
            ground.CONST.L_f = [];  % volumetric latent heat of fusion, freezing
            ground.CONST.c_w = [];  % volumetric heat capacity water
            ground.CONST.c_i = [];  % volumetric heat capacity ice
            ground.CONST.c_o = [];  % volumetric heat capacity organic
            ground.CONST.c_m = [];  % volumetric heat capacity mineral
            
            ground.CONST.k_a = [];   % thermal conductivity air
            ground.CONST.k_w = [];   % thermal conductivity water
            ground.CONST.k_i = [];   % thermal conductivity ice 
            ground.CONST.k_o = [];   % thermal conductivity organic 
            ground.CONST.k_m = [];   % thermal conductivity mineral 
            
            ground.CONST.sigma = []; %Stefan-Boltzmann constant
            ground.CONST.kappa = []; % von Karman constant
            ground.CONST.L_s = [];  %latent heat of sublimation, latent heat of evaporation handled in a dedicated function
            
            ground.CONST.cp = []; % specific heat capacity at constant pressure of air
            ground.CONST.g = []; % gravitational acceleration Earth surface
            
            ground.CONST.rho_w = [];   % water density
            ground.CONST.rho_i = [];   %ice density
        end

        function ground = convert_units(ground, tile)
            unit_converter = str2func(tile.PARA.unit_conversion_class);
            unit_converter = unit_converter();
            ground = convert_normal(unit_converter, ground, tile);
        end
        
        function ground = finalize_init(ground, tile)

            if isempty(ground.PARA.conductivity_function) || sum(isnan(ground.PARA.conductivity_function))>0
                ground.PARA.conductivity_function = 'conductivity_mixing_squares_heatPump';
            end
            
            total_area = ground.STATVAR.area(1,1);

            distance_between_circular_cells = (ground.PARA.distance_between_pipes-ground.PARA.pipe_diameter) ./ ground.PARA.number_of_circular_cells ./ 2;
            area_pipes = pi/4 .* ground.PARA.pipe_diameter.^2 .* ground.PARA.number_of_pipes;
            contact_length_pipes = pi .* ground.PARA.pipe_diameter .* ground.PARA.number_of_pipes;
            area_circular_cells = area_pipes;
            contact_length_circular_cells = contact_length_pipes;
            lateral_thickness = distance_between_circular_cells;
            for i=1:ground.PARA.number_of_circular_cells
                area_circular_cells = [pi/4.*(ground.PARA.pipe_diameter+2*i*distance_between_circular_cells).^2.*ground.PARA.number_of_pipes-area_circular_cells(:, 1) area_circular_cells];
                contact_length_circular_cells = [pi.*(ground.PARA.pipe_diameter+2*i*distance_between_circular_cells).*ground.PARA.number_of_pipes contact_length_circular_cells];      
                lateral_thickness = [lateral_thickness distance_between_circular_cells];
            end

            ground.STATVAR.area = repmat([total_area - sum(area_circular_cells) area_circular_cells], size(ground.STATVAR.layerThick,1),1); %first column: outside; last column: pipes
            ground.STATVAR.contact_length = repmat([contact_length_circular_cells], size(ground.STATVAR.layerThick,1),1); %first column: outside; last column: pipes
            ground.STATVAR.lateral_thickness = repmat([lateral_thickness lateral_thickness(:,end).*0], size(ground.STATVAR.layerThick,1),1);

            waterIce = [];
            organic = [];
            mineral = [];
            cooling_fluid=[];
            for i=1:size(ground.STATVAR.layerThick,1)
                if ~ground.STATVAR.heatPump_active(i,1)
                    waterIce = [waterIce; ground.STATVAR.waterIce(i,1).*ground.STATVAR.area(i,:)./total_area];
                    organic = [organic; ground.STATVAR.organic(i,1).*ground.STATVAR.area(i,:)./total_area];
                    mineral = [mineral; ground.STATVAR.mineral(i,1).*ground.STATVAR.area(i,:)./total_area];
                    cooling_fluid = [cooling_fluid; ground.STATVAR.area(i,:).*0];
                else
                    waterIce = [waterIce; [ground.STATVAR.waterIce(i,1).*ground.STATVAR.area(i,1:end-1)./total_area 0]]; %take out the part taken by the pipes
                    organic = [organic; [ground.STATVAR.organic(i,1).*ground.STATVAR.area(i,1:end-1)./total_area 0]];
                    mineral = [mineral; [ground.STATVAR.mineral(i,1).*ground.STATVAR.area(i,1:end-1)./total_area 0]];
                    cooling_fluid = [cooling_fluid; [ground.STATVAR.area(i,1:end-1).*0 ground.STATVAR.area(i,end).*ground.STATVAR.layerThick(i,end)]];

                end
            end
            ground.STATVAR.waterIce = waterIce; 
            ground.STATVAR.mineral = mineral;
            ground.STATVAR.organic = organic; 
            ground.STATVAR.cooling_fluid = cooling_fluid;
            ground.STATVAR.T = repmat(ground.STATVAR.T, 1, size(ground.STATVAR.waterIce,2));
            ground.STATVAR.layerThick = repmat(ground.STATVAR.layerThick, 1, size(ground.STATVAR.waterIce,2));

            ground = get_E_freeW_heatPump(ground);

            ground.TEMP.upper_heatPunp_cell = find(ground.STATVAR.heatPump_active(:,1)==1, 1, 'first');
            ground.TEMP.lower_heatPunp_cell = find(ground.STATVAR.heatPump_active(:,1)==1, 1, 'last');
            ground.STATVAR.incoming_fluid_T = ground.PARA.incoming_fluid_T;
            ground.STATVAR.fluid_speed = ground.PARA.fluid_speed;

            ground.STATVAR.Lstar = -100;
            ground.STATVAR.Qh = 0;
            ground.STATVAR.Qe = 0;
            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
            
        end
        
        function ground = finalize_init2(ground, tile)

            ground = get_E_freeW_heatPump(ground);
            
        end
        
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)
            forcing = tile.FORCING;
            T = ground.STATVAR.T;
            ground.STATVAR.T = sum(ground.STATVAR.T(1,:).*ground.STATVAR.area(1,:))./ sum(ground.STATVAR.area(1,:));
            ground = surface_energy_balance(ground, forcing);
            ground.STATVAR.T = T;
        end
        
        function [ground, S_up] = penetrate_SW(ground, S_down)  %mandatory function when used with class that features SW penetration 
            [ground, S_up] = penetrate_SW_no_transmission(ground, S_down);
        end
        
        function [ground, S_up] = penetrate_SW_PARENT(ground, S_down)  %mandatory function when used with class that features SW penetration
            [ground, S_up] = penetrate_SW_no_transmission(ground, S_down);
        end
        
        function ground = get_boundary_condition_l(ground, tile)
            forcing = tile.FORCING;
            ground.TEMP.F_lb = forcing.PARA.heatFlux_lb .* ground.STATVAR.area(end,:);
            ground.TEMP.d_energy(end,:) = ground.TEMP.d_energy(end,:) + ground.TEMP.F_lb;
        end
        
        function ground = get_derivatives_prognostic(ground, tile)

            ground = get_derivative_energy_heatPump(ground);

            %advection of cooling fluid
            ground.STATVAR.F_cooling = ground.TEMP.d_energy(:,end).*0;
            range = ground.TEMP.upper_heatPunp_cell:ground.TEMP.lower_heatPunp_cell;
            if ground.STATVAR.fluid_speed > 0
                ground.STATVAR.F_cooling(ground.TEMP.upper_heatPunp_cell,1) = ground.STATVAR.F_cooling(ground.TEMP.upper_heatPunp_cell,1) + ground.STATVAR.incoming_fluid_T .* ground.PARA.heatCapacity_fluid .* ground.STATVAR.fluid_speed;
                ground.STATVAR.F_cooling(range,1) = ground.STATVAR.F_cooling(range,1) - ground.STATVAR.T(range,end) .* ground.PARA.heatCapacity_fluid .* ground.STATVAR.fluid_speed;
                ground.STATVAR.F_cooling(range(2:end),1) = ground.STATVAR.F_cooling(range(2:end),1) + ground.STATVAR.T(range(1:end-1),end) .* ground.PARA.heatCapacity_fluid .* ground.STATVAR.fluid_speed;
                ground.STATVAR.F_cooling_net  = -(ground.STATVAR.incoming_fluid_T - ground.STATVAR.T(ground.TEMP.lower_heatPunp_cell,end)) .* ground.PARA.heatCapacity_fluid .* ground.STATVAR.fluid_speed .* ground.STATVAR.area(1,end);
                ground.STATVAR.outgoing_fluid_T = ground.STATVAR.T(ground.TEMP.lower_heatPunp_cell,end);
            else
                ground.STATVAR.F_cooling(ground.TEMP.lower_heatPunp_cell,1) = ground.STATVAR.F_cooling(ground.TEMP.lower_heatPunp_cell,1) - ground.STATVAR.incoming_fluid_T .* ground.PARA.heatCapacity_fluid .* ground.STATVAR.fluid_speed;
                ground.STATVAR.F_cooling(range,1) = ground.STATVAR.F_cooling(range,1) + ground.STATVAR.T(range,end) .* ground.PARA.heatCapacity_fluid .* ground.STATVAR.fluid_speed;
                ground.STATVAR.F_cooling(range(1:end-1),1) = ground.STATVAR.F_cooling(range(1:end-1),1) - ground.STATVAR.T(range(2:end),end) .* ground.PARA.heatCapacity_fluid .* ground.STATVAR.fluid_speed;
                ground.STATVAR.F_cooling_net  = (ground.STATVAR.incoming_fluid_T - ground.STATVAR.T(ground.TEMP.upper_heatPunp_cell,end)) .* ground.PARA.heatCapacity_fluid .* ground.STATVAR.fluid_speed .* ground.STATVAR.area(1,end);
                ground.STATVAR.outgoing_fluid_T = ground.STATVAR.T(ground.TEMP.upper_heatPunp_cell,end);
            end
            
            ground.STATVAR.F_cooling = ground.STATVAR.F_cooling .* ground.STATVAR.area(:,end) ; %flow speed is flow speed in each indiuvidual pipe, but area is for all pipes
            ground.TEMP.d_energy(:,end) = ground.TEMP.d_energy(:,end) + ground.STATVAR.F_cooling;
           
        end
        
        function timestep = get_timestep(ground, tile) 
           timestep = min(get_timestep_heat_coduction(ground));
        end
        
        function ground = advance_prognostic(ground, tile)           
            timestep = tile.timestep;
            %energy
            ground.STATVAR.energy = ground.STATVAR.energy + timestep .* ground.TEMP.d_energy;
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
            forcing = tile.FORCING;
            ground = L_star(ground, forcing);
        end
       
        function ground = compute_diagnostic(ground, tile)
            ground = get_T_water_freeW_heatPump(ground);
            ground = conductivity(ground);
            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
        end
        
        function ground = check_trigger(ground, tile)
           %do nothing 
        end
        

        
        
        %-----non-mandatory functions-------
        function ground = surface_energy_balance(ground, forcing)
            ground.STATVAR.Lout = (1-ground.PARA.epsilon) .* forcing.TEMP.Lin + ground.PARA.epsilon .* ground.CONST.sigma .* (ground.STATVAR.T + 273.15).^4;
            ground.STATVAR.Sout = ground.PARA.albedo .*  forcing.TEMP.Sin;
            ground.STATVAR.Qh = Q_h(ground, forcing);
            ground.STATVAR.Qe = Q_eq(ground, forcing);
            
            ground.TEMP.F_ub = (forcing.TEMP.Sin + forcing.TEMP.Lin - ground.STATVAR.Lout - ground.STATVAR.Sout - ground.STATVAR.Qh - ground.STATVAR.Qe) .* ground.STATVAR.area(1,:);
            ground.TEMP.d_energy(1,:) = ground.TEMP.d_energy(1,:) + ground.TEMP.F_ub;
        end
        
        function ground = conductivity(ground)
            conductivity_function = str2func(ground.PARA.conductivity_function);
            ground = conductivity_function(ground);
        end

        function gridcell_variables = get_gridcell_variables(ground)
            gridcell_variables = {'waterIce'; 'mineral'; 'organic'; 'T'; 'energy'; 'area'; 'layerThick';  'air'; 'water'; 'ice'};

        end

        %-------------restore_from_out-----
        function ground = reset_from_OUT(ground, tile)
            ground.TEMP.d_energy = ground.TEMP.d_energy.*0; 
        end
        
        
        %-----LATERAL-------------------
        
        function ground = lateral_push_groundHeatPump(ground, lateral)

            if lateral.PARA.target_T_inside - lateral.STATVAR.Tair > 5 %heating only if Tair is 5 degreeC colder than air
                heating_power_need_house = (lateral.PARA.target_T_inside - lateral.STATVAR.Tair) ./ lateral.TEMP.R_house;
                ground.STATVAR.fluid_speed = ground.PARA.fluid_speed; %assumed constant
            
                T_vaporization = ground.STATVAR.outgoing_fluid_T - 5; %fluid coming from ground must be 5 degrees warmer than the vaporization T of the heat pump coolant, taken from internet (10 degreeC for air heat pumps)
                current_COP = interp1(lateral.PARA.COP(:,1), lateral.PARA.COP(:,2), T_vaporization);
                power_extracted_from_ground = heating_power_need_house .* (1-1/current_COP);
                ground.STATVAR.incoming_fluid_T = ground.STATVAR.outgoing_fluid_T - power_extracted_from_ground./ (ground.PARA.heatCapacity_fluid .* ground.STATVAR.fluid_speed .* ground.STATVAR.area(1,end)); %this is feedback to the ground class
                
                lateral.STATVAR.heating_power_need_house = heating_power_need_house;
                lateral.STATVAR.T_vaporization_ground = T_vaporization;
                lateral.STATVAR.COP_ground = current_COP;
                lateral.STATVAR.electric_power_ground = heating_power_need_house - power_extracted_from_ground;
                lateral.STATVAR.power_extracted_from_ground = power_extracted_from_ground;
                lateral.STATVAR.incoming_fluid_T = ground.STATVAR.incoming_fluid_T;
                lateral.STATVAR.outgoing_fluid_T = ground.STATVAR.outgoing_fluid_T;
                lateral.STATVAR.running_cost_ground = lateral.STATVAR.electric_power_ground .* lateral.STATVAR.electricity_price;

                %air heat Pump for comparison
                T_vaporization = lateral.STATVAR.Tair - 10; %fluid coming from ground must be 5 degrees warmer than the vaporization T of the heat pump coolant, taken from internet (10 degreeC for air heat pumps)
                current_COP = interp1(lateral.PARA.COP(:,1), lateral.PARA.COP(:,2), T_vaporization);
                power_extracted_from_air = heating_power_need_house .* (1-1/current_COP);
                lateral.STATVAR.T_vaporization_air = T_vaporization;
                lateral.STATVAR.COP_air = current_COP;
                lateral.STATVAR.electric_power_air = heating_power_need_house - power_extracted_from_air;
                lateral.STATVAR.power_extracted_from_air = power_extracted_from_air;
                lateral.STATVAR.running_cost_air = lateral.STATVAR.electric_power_air .* lateral.STATVAR.electricity_price;

            else
                ground.STATVAR.fluid_speed = 0;

                lateral.STATVAR.heating_power_need_house = 0;
                lateral.STATVAR.T_vaporization_ground = NaN;
                lateral.STATVAR.COP_ground = NaN;
                lateral.STATVAR.electric_power_ground = NaN;
                lateral.STATVAR.power_extracted_from_ground = NaN;
                lateral.STATVAR.incoming_fluid_T = NaN;
                lateral.STATVAR.outgoing_fluid_T = NaN;

                lateral.STATVAR.T_vaporization_air = NaN;
                lateral.STATVAR.COP_air = NaN;
                lateral.STATVAR.electric_power = NaN;
                lateral.STATVAR.power_extracted_from_air = NaN;
                lateral.STATVAR.running_cost_ground = 0;
                lateral.STATVAR.running_cost_air = 0;
            end

        end

        function ground = lateral_push_groundHeatPump_recharge(ground, lateral)

            
            if lateral.PARA.target_T_inside - lateral.STATVAR.Tair > 5 %heating only if Tair is 5 degreeC colder than air
                heating_power_need_house = (lateral.PARA.target_T_inside - lateral.STATVAR.Tair) ./ lateral.TEMP.R_house;
                lateral.STATVAR.heating_power_need_house = heating_power_need_house;

                range = ground.TEMP.upper_heatPunp_cell:ground.TEMP.lower_heatPunp_cell;
                if  lateral.STATVAR.electricity_price < 0.5*lateral.TEMP.average_electricity_price && mean(ground.STATVAR.T(range,end)) < 10 %cheap electricity and ground is cold -> heating house and recharging ground from air

                        T_vaporization = lateral.STATVAR.Tair - 10; %fluid coming from ground must be 5 degrees warmer than the vaporization T of the heat pump coolant, taken from internet (10 degreeC for air heat pumps)
                        current_COP = interp1(lateral.PARA.COP(:,1), lateral.PARA.COP(:,2), T_vaporization);
                        
                        ground.STATVAR.incoming_fluid_T = 35; %same T as for house
                        
                        ground.STATVAR.fluid_speed = min(ground.PARA.fluid_speed,  (ground.PARA.max_power .* current_COP - heating_power_need_house) ./ ...
                            ((ground.STATVAR.incoming_fluid_T - ground.STATVAR.outgoing_fluid_T) .* ground.PARA.heatCapacity_fluid .* ground.STATVAR.area(1,end)));
                        power_injected_in_ground = (ground.STATVAR.incoming_fluid_T - ground.STATVAR.outgoing_fluid_T) .* ground.PARA.heatCapacity_fluid .* ground.STATVAR.fluid_speed .* ground.STATVAR.area(1,end); 

                        ground.STATVAR.fluid_speed = -ground.STATVAR.fluid_speed; %reversing flow, so that strongest heating is in the bottom -> could down-regulate to hold max-power of heatpump

                        power_extracted_from_air = (heating_power_need_house+power_injected_in_ground) .* (1-1/current_COP);

                        lateral.STATVAR.T_vaporization = T_vaporization;
                        lateral.STATVAR.COP = current_COP;
                        lateral.STATVAR.electric_power = heating_power_need_house + power_injected_in_ground - power_extracted_from_air;
                        lateral.STATVAR.power_extracted = power_extracted_from_air;
                        lateral.STATVAR.incoming_fluid_T = ground.STATVAR.incoming_fluid_T;
                        lateral.STATVAR.outgoing_fluid_T = ground.STATVAR.outgoing_fluid_T;
                        lateral.STATVAR.running_cost = lateral.STATVAR.electric_power .* lateral.STATVAR.electricity_price;
                        lateral.STATVAR.run_state = 10;
                else %check whether air or ground heat pump is more efficient
                    if lateral.STATVAR.Tair - 10 > ground.STATVAR.outgoing_fluid_T - 5 %use air heat pump
                        
                        ground.STATVAR.fluid_speed = 0;
                        %air heat Pump for comparison
                        T_vaporization = lateral.STATVAR.Tair - 10; %fluid coming from ground must be 5 degrees warmer than the vaporization T of the heat pump coolant, taken from internet (10 degreeC for air heat pumps)
                        current_COP = interp1(lateral.PARA.COP(:,1), lateral.PARA.COP(:,2), T_vaporization);
                        power_extracted_from_air = heating_power_need_house .* (1-1/current_COP);
                        lateral.STATVAR.T_vaporization = T_vaporization;
                        lateral.STATVAR.COP = current_COP;
                        lateral.STATVAR.electric_power = heating_power_need_house - power_extracted_from_air;
                        lateral.STATVAR.power_extracted = power_extracted_from_air;
                        lateral.STATVAR.incoming_fluid_T = NaN;
                        lateral.STATVAR.outgoing_fluid_T = NaN;
                        lateral.STATVAR.running_cost = lateral.STATVAR.electric_power .* lateral.STATVAR.electricity_price;
                        lateral.STATVAR.run_state = 1;
                    else %use ground heat pump
                        ground.STATVAR.fluid_speed = ground.PARA.fluid_speed; %assumed constant

                        T_vaporization = ground.STATVAR.outgoing_fluid_T - 5; %fluid coming from ground must be 5 degrees warmer than the vaporization T of the heat pump coolant, taken from internet (10 degreeC for air heat pumps)
                        current_COP = interp1(lateral.PARA.COP(:,1), lateral.PARA.COP(:,2), T_vaporization);
                        power_extracted_from_ground = heating_power_need_house .* (1-1/current_COP);
                        ground.STATVAR.incoming_fluid_T = ground.STATVAR.outgoing_fluid_T - power_extracted_from_ground./ (ground.PARA.heatCapacity_fluid .* ground.STATVAR.fluid_speed .* ground.STATVAR.area(1,end)); %this is feedback to the ground class

                        lateral.STATVAR.T_vaporization_ground = T_vaporization;
                        lateral.STATVAR.COP = current_COP;
                        lateral.STATVAR.electric_power = heating_power_need_house - power_extracted_from_ground;
                        lateral.STATVAR.power_extracted = power_extracted_from_ground;
                        lateral.STATVAR.incoming_fluid_T = ground.STATVAR.incoming_fluid_T;
                        lateral.STATVAR.outgoing_fluid_T = ground.STATVAR.outgoing_fluid_T;
                        lateral.STATVAR.running_cost = lateral.STATVAR.electric_power .* lateral.STATVAR.electricity_price;
                        lateral.STATVAR.run_state = 2;
                    end

                end

            else
                heating_power_need_house = 0;
                ground.STATVAR.fluid_speed = 0;

                lateral.STATVAR.heating_power_need_house = 0;
                lateral.STATVAR.T_vaporization = NaN;
                lateral.STATVAR.COP = NaN;
                lateral.STATVAR.electric_power = NaN;
                lateral.STATVAR.power_extracted = NaN;
                lateral.STATVAR.incoming_fluid_T = NaN;
                lateral.STATVAR.outgoing_fluid_T = NaN;
                lateral.STATVAR.running_cost = 0;
                lateral.STATVAR.run_state = 0;
            end
        end

        function ground = lateral_push_groundHeatPump_recharge2(ground, lateral)

            
            if lateral.PARA.target_T_inside - lateral.STATVAR.Tair > 5 %heating only if Tair is 5 degreeC colder than air
                heating_power_need_house = (lateral.PARA.target_T_inside - lateral.STATVAR.Tair) ./ lateral.TEMP.R_house;
                lateral.STATVAR.heating_power_need_house = heating_power_need_house;

                range = ground.TEMP.upper_heatPunp_cell:ground.TEMP.lower_heatPunp_cell;
                if  lateral.STATVAR.electricity_price < 0.5*lateral.TEMP.average_electricity_price && mean(ground.STATVAR.T(range,end)) < 10 %cheap electricity and ground is cold -> heating house and recharging ground from air

                        T_vaporization = lateral.STATVAR.Tair - 10; %fluid coming from ground must be 5 degrees warmer than the vaporization T of the heat pump coolant, taken from internet (10 degreeC for air heat pumps)
                        current_COP_house = interp1(lateral.PARA.COP(:,1), lateral.PARA.COP(:,2), T_vaporization);
                        power_extracted_from_air = heating_power_need_house .* (1-1/current_COP_house);
                        lateral.STATVAR.electric_power = heating_power_need_house - power_extracted_from_air;

                        max_power_injected_in_ground = max(0, ground.PARA.max_power - (heating_power_need_house - power_extracted_from_air));

                        current_COP_ground = interp1(lateral.PARA.COP2(:,1), lateral.PARA.COP2(:,2), T_vaporization);
                        ground.STATVAR.incoming_fluid_T = 15; %lower T than for house to increase efficiency
                        ground.STATVAR.fluid_speed = min(ground.PARA.fluid_speed,  max_power_injected_in_ground ./ ...
                            ((ground.STATVAR.incoming_fluid_T - ground.STATVAR.outgoing_fluid_T) .* ground.PARA.heatCapacity_fluid .* ground.STATVAR.area(1,end)));
                        power_injected_in_ground = (ground.STATVAR.incoming_fluid_T - ground.STATVAR.outgoing_fluid_T) .* ground.PARA.heatCapacity_fluid .* ground.STATVAR.fluid_speed .* ground.STATVAR.area(1,end); 

                        ground.STATVAR.fluid_speed = -ground.STATVAR.fluid_speed; %reversing flow, so that strongest heating is in the bottom -> could down-regulate to hold max-power of heatpump
                        power_extracted_from_air2 = power_injected_in_ground .* (1-1/current_COP_ground);
                        lateral.STATVAR.electric_power = lateral.STATVAR.electric_power + (power_injected_in_ground - power_extracted_from_air2);

                        lateral.STATVAR.T_vaporization = T_vaporization;
                        lateral.STATVAR.COP = current_COP_house;
                        lateral.STATVAR.COP2 = current_COP_ground;
                        lateral.STATVAR.power_extracted = power_extracted_from_air+power_extracted_from_air2;
                        lateral.STATVAR.incoming_fluid_T = ground.STATVAR.incoming_fluid_T;
                        lateral.STATVAR.outgoing_fluid_T = ground.STATVAR.outgoing_fluid_T;
                        lateral.STATVAR.running_cost = lateral.STATVAR.electric_power .* lateral.STATVAR.electricity_price;
                        lateral.STATVAR.run_state = 10;
                else %check whether air or ground heat pump is more efficient
                    if lateral.STATVAR.Tair - 10 > ground.STATVAR.outgoing_fluid_T - 5 %use air heat pump
                        
                        ground.STATVAR.fluid_speed = 0;
                        %air heat Pump for comparison
                        T_vaporization = lateral.STATVAR.Tair - 10; %fluid coming from ground must be 5 degrees warmer than the vaporization T of the heat pump coolant, taken from internet (10 degreeC for air heat pumps)
                        current_COP = interp1(lateral.PARA.COP(:,1), lateral.PARA.COP(:,2), T_vaporization);
                        power_extracted_from_air = heating_power_need_house .* (1-1/current_COP);
                        lateral.STATVAR.T_vaporization = T_vaporization;
                        lateral.STATVAR.COP = current_COP;
                        lateral.STATVAR.COP2 = NaN;
                        lateral.STATVAR.electric_power = heating_power_need_house - power_extracted_from_air;
                        lateral.STATVAR.power_extracted = power_extracted_from_air;
                        lateral.STATVAR.incoming_fluid_T = NaN;
                        lateral.STATVAR.outgoing_fluid_T = NaN;
                        lateral.STATVAR.running_cost = lateral.STATVAR.electric_power .* lateral.STATVAR.electricity_price;
                        lateral.STATVAR.run_state = 1;
                    else %use ground heat pump
                        ground.STATVAR.fluid_speed = ground.PARA.fluid_speed; %assumed constant

                        T_vaporization = ground.STATVAR.outgoing_fluid_T - 5; %fluid coming from ground must be 5 degrees warmer than the vaporization T of the heat pump coolant, taken from internet (10 degreeC for air heat pumps)
                        current_COP = interp1(lateral.PARA.COP(:,1), lateral.PARA.COP(:,2), T_vaporization);
                        power_extracted_from_ground = heating_power_need_house .* (1-1/current_COP);
                        ground.STATVAR.incoming_fluid_T = ground.STATVAR.outgoing_fluid_T - power_extracted_from_ground./ (ground.PARA.heatCapacity_fluid .* ground.STATVAR.fluid_speed .* ground.STATVAR.area(1,end)); %this is feedback to the ground class

                        lateral.STATVAR.T_vaporization_ground = T_vaporization;
                        lateral.STATVAR.COP = current_COP;
                        lateral.STATVAR.COP2 = NaN;
                        lateral.STATVAR.electric_power = heating_power_need_house - power_extracted_from_ground;
                        lateral.STATVAR.power_extracted = power_extracted_from_ground;
                        lateral.STATVAR.incoming_fluid_T = ground.STATVAR.incoming_fluid_T;
                        lateral.STATVAR.outgoing_fluid_T = ground.STATVAR.outgoing_fluid_T;
                        lateral.STATVAR.running_cost = lateral.STATVAR.electric_power .* lateral.STATVAR.electricity_price;
                        lateral.STATVAR.run_state = 2;
                    end

                end

            else
                heating_power_need_house = 0;
                ground.STATVAR.fluid_speed = 0;

                lateral.STATVAR.heating_power_need_house = 0;
                lateral.STATVAR.T_vaporization = NaN;
                lateral.STATVAR.COP = NaN;
                lateral.STATVAR.COP2 = NaN;
                lateral.STATVAR.electric_power = NaN;
                lateral.STATVAR.power_extracted = NaN;
                lateral.STATVAR.incoming_fluid_T = NaN;
                lateral.STATVAR.outgoing_fluid_T = NaN;
                lateral.STATVAR.running_cost = 0;
                lateral.STATVAR.run_state = 0;
            end
        end




        %-------LAT3D_HEAT-------------
        function ground = lateral3D_pull_heat(ground, lateral)
            ground = lateral3D_pull_heat_simple(ground, lateral);
        end
        
        function ground = lateral3D_push_heat(ground, lateral)
            ground = lateral3D_push_heat_simple(ground, lateral);
        end
        
        %----ADJUST_STRATIGRAPHY-------------
        %needed to change the stratigraphy during data assimilation
        
        function ground = adjust_stratigraphy(ground, tile)
            ground = adjust_stratigraphy_waterIce_fixed(ground, tile);
        end
        
                 
         %-------------param file generation-----
         function ground = param_file_info(ground)
             ground = param_file_info@BASE(ground);
             %ground = provide_PARA(ground);
             
             ground.PARA.class_category = 'GROUND';
             
             %ground.PARA.options = [];
             ground.PARA.STATVAR = {'waterIce' 'mineral' 'organic' 'T'};
             
             ground.PARA.default_value.albedo = {0.2};
             ground.PARA.comment.albedo = {'surface albedo [-]'};
             
             ground.PARA.default_value.epsilon = {0.99};
             ground.PARA.comment.epsilon = {'surface emissivity [-]'};
             
             ground.PARA.default_value.z0 = {0.01};
             ground.PARA.comment.z0 = {'roughness length [m]'};
             
             ground.PARA.default_value.rs = {0};
             ground.PARA.comment.rs ={'surface resistance against evapotranspiration [sec/m]'};
             
             ground.PARA.default_value.conductivity_function = {'conductivity_mixing_squares'};
             ground.PARA.comment.conductivity_function = {'function employed to calculate thermal conductivity, leave empty for default'};
             
             ground.PARA.default_value.dt_max = {3600};
             ground.PARA.comment.dt_max = {'maximum possible timestep [sec]'};
             
             ground.PARA.default_value.dE_max = {50000};
             ground.PARA.comment.dE_max = {'maximum possible energy change per timestep [J/m3]'};
        end
        
        %----inherited Tier 1 functions ------------
        
%         function ground = get_derivative_energy(ground)
%            ground = get_derivative_energy@HEAT_CONDUCTION(ground); 
%         end
%         
%         function ground = conductivity_mixing_squares(ground)
%             ground = conductivity_mixing_squares@HEAT_CONDUCTION(ground);
%         end
%         
%         function flux = Q_h(ground, forcing)
%            flux = Q_h@SEB(ground, forcing);
%         end
%     
%         function flux = Q_eq(ground, forcing)
%             flux = Q_eq@SEB(ground, forcing);
%         end
%         
%         function timestep = get_timestep_heat_coduction(ground)
%             timestep = get_timestep_heat_coduction@HEAT_CONDUCTION(ground);
%         end
%         
%         function ground = L_star(ground, forcing)
%            ground = L_star@SEB(ground, forcing); 
%         end
%         
%         function [ground, S_up] = penetrate_SW_no_transmission(ground, S_down)
%             [ground, S_up] = penetrate_SW_no_transmission@SEB(ground, S_down);
%         end
%         
%         function ground = get_T_water_freeW(ground)
%             ground = get_T_water_freeW@HEAT_CONDUCTION(ground);
%         end
    end
    
end
