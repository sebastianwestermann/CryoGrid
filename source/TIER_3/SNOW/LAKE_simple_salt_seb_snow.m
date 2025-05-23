%========================================================================
% CryoGrid GROUND class LAKE_simple_salt_seb_snow
% static water body with heat conduction, free water freeze curve, surface
% energy balance
% representation of frozen water body, works in concert with
% LAKE_simple_unfrozen_seb for unfrozen water body
% S. Westermann, October 2020
%========================================================================

classdef LAKE_simple_salt_seb_snow < LAKE_simple_salt_seb
    
    properties
        CHILD
        IA_CHILD
    end
    
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        function ground = initialize_from_LAKE_previous_season(ground, LAKE_simple_unfrozen)
            ground = initialize_from_LAKE_previous_season@LAKE_simple_salt_seb(ground, LAKE_simple_unfrozen);
            ground.CHILD = 0; % no snow
            ground.IA_CHILD = 0;
        end
        
        function ground = provide_PARA(ground)
            ground = provide_PARA@LAKE_simple_salt_seb(ground);
        end
        
        function ground = provide_CONST(ground)
            ground = provide_CONST@LAKE_simple_salt_seb(ground);
        end
        
        function ground = provide_STATVAR(ground)
            ground = provide_STATVAR@LAKE_simple_salt_seb(ground);
        end
        
        function ground = finalize_init(ground, tile)
            ground = finalize_init@LAKE_simple_salt_seb(ground, tile);
            ground.CHILD = 0; % no snow
            ground.IA_CHILD = 0;
        end
        
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)
            forcing = tile.FORCING;
            if ground.CHILD == 0  %CHILD does not exist
                ground = get_boundary_condition_u@LAKE_simple_salt_seb(ground, tile); %call the native function for the ground class
                
                if forcing.TEMP.snowfall > 0  %create CHILD
%                     CURRENT = ground.PREVIOUS;  %go to Top() and get the stored SNOW class
%                     while ~strcmp(class(CURRENT), 'Top')
%                         CURRENT = CURRENT.PREVIOUS;
%                     end
%                     ground.CHILD = copy(CURRENT.STORE.SNOW);
                    ground.CHILD = copy(tile.STORE.SNOW);
                    ground.CHILD.PARENT = ground;
                    ground.CHILD.NEXT = ground;
                    ground.IA_CHILD = get_IA_class(class(ground.CHILD), class(ground));
                    ground.IA_CHILD.PREVIOUS = ground.CHILD;
                    ground.IA_CHILD.NEXT = ground; %SNOW CHILD created
                    
                    ground.CHILD = get_boundary_condition_u_create_CHILD(ground.CHILD, tile);  %initialize with fresh snowfall
                end
            else %CHILD exists
                total_area = ground.STATVAR.area(1,1); %store the total area of the ground
                ground.STATVAR.area(1,1) = ground.STATVAR.area(1,1) - ground.CHILD.STATVAR.area(1,1); %replace by snow-free area
                
                ground.CHILD.STATVAR.Lstar = ground.STATVAR.Lstar;
                
                ground.CHILD = get_boundary_condition_u_CHILD(ground.CHILD, tile);
                ground = get_boundary_condition_u@LAKE_simple_salt_seb(ground, tile);
                
                get_IA_CHILD_boundary_condition_u(ground.IA_CHILD, tile);
                %call designated mandatory function for CHILD-PARENT interactions in
                %the IA class governing IA between SNOW and GROUND
                
                ground.STATVAR.Lout = (ground.STATVAR.area(1,1) .* ground.STATVAR.Lout + ground.CHILD.STATVAR.area .* ground.CHILD.STATVAR.Lout) ./ total_area; %mix the surface heat fluxes from snow and ground
                ground.STATVAR.Sout = (ground.STATVAR.area(1,1) .* ground.STATVAR.Sout + ground.CHILD.STATVAR.area .* ground.CHILD.STATVAR.Sout) ./ total_area;
                ground.STATVAR.Qh = (ground.STATVAR.area(1,1) .* ground.STATVAR.Qh + ground.CHILD.STATVAR.area .* ground.CHILD.STATVAR.Qh) ./ total_area;
                ground.STATVAR.Qe = (ground.STATVAR.area(1,1) .* ground.STATVAR.Qe + ground.CHILD.STATVAR.area .* ground.CHILD.STATVAR.Qe) ./ total_area;
                
                ground.STATVAR.area(1,1) = total_area; %reassign the true area of ground
            end
        end
        
        function ground = get_boundary_condition_l(ground, tile)
            ground = get_boundary_condition_l@LAKE_simple_salt_seb(ground, tile);
        end
        
        function [ground, S_up] = penetrate_SW(ground, S_down)  %mandatory function when used with class that features SW penetration
            [ground, S_up] = penetrate_SW@LAKE_simple_salt_seb(ground, S_down);
        end

        function [ground, S_up] = penetrate_SW_PARENT(ground, S_down)  %mandatory function when used with class that features SW penetration
            [ground, S_up] = penetrate_SW_PARENT@LAKE_simple_salt_seb(ground, S_down);
        end
        
        function ground = get_derivatives_prognostic(ground, tile)
            if ground.CHILD == 0
                ground = get_derivatives_prognostic@LAKE_simple_salt_seb(ground, tile); %call normal function
            else
                ground.CHILD = get_derivatives_prognostic_CHILD(ground.CHILD, tile);
                ground = get_derivatives_prognostic@LAKE_simple_salt_seb(ground, tile);
            end
        end
        
        function timestep = get_timestep(ground, tile)
            if ground.CHILD == 0
                timestep = get_timestep@LAKE_simple_salt_seb(ground, tile);
            else
                timestep_snow = get_timestep_CHILD(ground.CHILD, tile);
                timestep_ground =  get_timestep@LAKE_simple_salt_seb(ground, tile);
                timestep = timestep_ground + double(timestep_snow > 0 && timestep_snow < timestep_ground) .* (timestep_snow - timestep_ground);
            end
        end
        
        function ground = advance_prognostic(ground, tile)
            if ground.CHILD == 0
                ground =  advance_prognostic@LAKE_simple_salt_seb(ground, tile);
            else
                ground.CHILD = advance_prognostic_CHILD(ground.CHILD, tile);
                ground =  advance_prognostic@LAKE_simple_salt_seb(ground, tile);
            end
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
            forcing = tile.FORCING;
            ground = L_star(ground, forcing);
        end
        
        function ground = compute_diagnostic(ground, tile)
            if ground.CHILD == 0
                ground = compute_diagnostic@LAKE_simple_salt_seb(ground, tile);
            else
                ground = compute_diagnostic@LAKE_simple_salt_seb(ground, tile);
                ground.CHILD = compute_diagnostic_CHILD(ground.CHILD, tile);
                
            end
        end
        
        function ground = check_trigger(ground, tile)
            
            %snow trigger
            if ground.CHILD ~= 0
                %delete CHILD
                if ground.CHILD.STATVAR.area ./ ground.STATVAR.area(1,1) < 1e-6 %cutoff to get rid of remaining snow
                    ground.CHILD = 0;
                    ground.IA_CHILD = 0;
                    %make SNOW CHILD full class
                elseif ground.CHILD.STATVAR.area ./ ground.STATVAR.area(1,1) > 1
                    %transforms dimensions and STAVAR
                    snow_volume = ground.CHILD.STATVAR.area .* ground.CHILD.STATVAR.layerThick;
                    ground.CHILD.STATVAR.area = ground.STATVAR.area(1,1);
                    ground.CHILD.STATVAR.layerThick = snow_volume ./ ground.CHILD.STATVAR.area;
                    %ground.CHILD = compute_diagnostic(ground.CHILD, forcing); %splits snow in 2 grid cells
                    
%                     %make snow a real class
%                     ground.CHILD.PARENT = 0;
%                     ground.CHILD.PREVIOUS = ground.PREVIOUS;
%                     ground.CHILD.NEXT = ground;
%                     ground.PREVIOUS.NEXT = ground.CHILD;
%                     ground.PREVIOUS = ground.CHILD;
%                     ground.CHILD = 0;
%                     ground.IA_PREVIOUS = ground.IA_CHILD;
%                     ground.PREVIOUS.IA_NEXT = ground.IA_CHILD;
%                     ground.IA_CHILD = 0;
                    
                    %make snow a real class
                    ground.CHILD.PARENT = 0;
                    ground.CHILD.PREVIOUS = ground.PREVIOUS;
                    ground.CHILD.NEXT = ground;
                    ground.PREVIOUS.NEXT = ground.CHILD;
                    ia_class = get_IA_class(class(ground.PREVIOUS), class(ground.CHILD));
                    ground.PREVIOUS.IA_NEXT = ia_class;
                    ground.CHILD.IA_PREVIOUS = ia_class;
                    ground.CHILD.IA_PREVIOUS.NEXT = ground.CHILD;
                    ground.CHILD.IA_PREVIOUS.PREVIOUS = ground.PREVIOUS;
                    finalize_init(ground.CHILD.IA_PREVIOUS, tile);
                    
                    ground.PREVIOUS = ground.CHILD;
                    ground.CHILD = 0;
                    ground.IA_PREVIOUS = ground.IA_CHILD; 
                    ground.PREVIOUS.IA_NEXT = ground.IA_CHILD;
                    ground.IA_CHILD = 0;
                end
            end
            
            %lake trigger
            dummy = check_trigger@LAKE_simple_salt_seb(ground, tile);
        end
        
                
        %----------
        %reset timestamp when changing TILES
        function ground = reset_timestamps(ground, tile)
            if ground.CHILD ~= 0
                ground.CHILD = reset_timestamps(ground.CHILD, tile);
            end
        end
        
    end
end