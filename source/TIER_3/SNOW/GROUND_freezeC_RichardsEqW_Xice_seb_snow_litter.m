%========================================================================
% CryoGrid GROUND class GROUND_freezeC_RichardsEqW_Xice_seb_snow_litter
% like GROUND_freezeC_RichardsEqW_Xice_seb_snow, but includes a layer of
% litter that emerges in fall and decreases until end of year
% R. B. Zweigel, February 2023
%========================================================================

classdef GROUND_freezeC_RichardsEqW_Xice_seb_snow_litter < GROUND_freezeC_RichardsEqW_Xice_seb_snow
    properties
        
    end
    
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        
        function ground = provide_PARA(ground)
            ground = provide_PARA@GROUND_freezeC_RichardsEqW_Xice_seb_snow(ground);
        end
        
        function ground = provide_CONST(ground)
            ground = provide_CONST@GROUND_freezeC_RichardsEqW_Xice_seb_snow(ground);
        end
        
        function ground = provide_STATVAR(ground)
            ground = provide_STATVAR@GROUND_freezeC_RichardsEqW_Xice_seb_snow(ground);
        end
        
        function ground = finalize_init(ground, tile)
            ground = finalize_init@GROUND_freezeC_RichardsEqW_Xice_seb_snow(ground, tile);
            ground.TEMP.litter = 0;
            ground.TEMP.target_layerThick = ground.STATVAR.layerThick(1); % Assumes we start without XwaterIce
        end
        
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)
            ground = get_boundary_condition_u@GROUND_freezeC_RichardsEqW_Xice_seb_snow(ground, tile);
        end
        
        
        function [ground, S_up] = penetrate_SW(ground, S_down)  %mandatory function when used with class that features SW penetration
            [ground, S_up] = penetrate_SW@GROUND_freezeC_RichardsEqW_Xice_seb_snow(ground, S_down);
        end
        
        function [ground, S_up] = penetrate_SW_PARENT(ground, S_down)  %mandatory function when used with class that features SW penetration
            [ground, S_up] = penetrate_SW_PARENT@GROUND_freezeC_RichardsEqW_Xice_seb_snow(ground, S_down);
        end
        
        function [ground, L_up] = penetrate_LW(ground, L_down)  %mandatory function when used with class that features SW penetration
            [ground, L_up] = penetrate_LW@GROUND_freezeC_RichardsEqW_Xice_seb_snow(ground, L_down);
        end
        
        
        function ground = get_boundary_condition_l(ground, tile)
            ground = get_boundary_condition_l@GROUND_freezeC_RichardsEqW_Xice_seb_snow(ground, tile);
        end
        
        
        function ground = get_derivatives_prognostic(ground, tile)
            ground = get_derivatives_prognostic@GROUND_freezeC_RichardsEqW_Xice_seb_snow(ground, tile);
        end
        
        function timestep = get_timestep(ground, tile)
            timestep = get_timestep@GROUND_freezeC_RichardsEqW_Xice_seb_snow(ground, tile);
        end
        
        function ground = advance_prognostic(ground, tile)
            ground =  advance_prognostic@GROUND_freezeC_RichardsEqW_Xice_seb_snow(ground, tile);
            
%             aboveClass = class(ground.PREVIOUS);
%             if ground.STATVAR.sublimation<0 && tile.timestep>0 && ~strcmp(aboveClass(1:4),'SNOW')
%                 potential_ice_volume = ground.STATVAR.layerThick(1).*ground.STATVAR.area(1) - ground.STATVAR.mineral(1) - ground.STATVAR.organic(1) - ground.STATVAR.XwaterIce(1);
%                 ice_available = ground.STATVAR.Xice(1) + min(ground.STATVAR.ice(1),.9.*potential_ice_volume);
%                 sublimation = max(ground.STATVAR.sublimation.*tile.timestep, -ice_available);
%                 
%                 ground.TEMP.sublimation_energy = ground.TEMP.sublimation_energy .* (sublimation./tile.timestep)/ground.STATVAR.sublimation;
%                 ground.STATVAR.sublimation = sublimation./tile.timestep;
%                 
%                 Xice_removed = max(sublimation,-ground.STATVAR.Xice(1));
%                 ground.STATVAR.XwaterIce(1) = max(0, ground.STATVAR.XwaterIce(1) + Xice_removed);
%                 ground.STATVAR.layerThick(1) = ground.STATVAR.layerThick(1) + Xice_removed./ground.STATVAR.area(1);
%                 ground.STATVAR.Xice(1) = max(0, ground.STATVAR.Xice(1) + Xice_removed);
%                 
%                 ice_removed = sublimation - Xice_removed;
%                 ground.STATVAR.waterIce(1) = ground.STATVAR.waterIce(1) + ice_removed;
%                 ground.STATVAR.ice(1) = ground.STATVAR.ice(1) + ice_removed;
%                 
%                 ground.STATVAR.energy(1) = ground.STATVAR.energy(1) + ground.TEMP.sublimation_energy.*tile.timestep;
%             end
            
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
            forcing = tile.FORCING;
            ground = L_star(ground, forcing);
        end
        
        function ground = compute_diagnostic(ground, tile)
%             ground = compute_diagnostic@GROUND_freezeC_RichardsEqW_Xice_seb_snow(ground, tile);
            if ground.CHILD == 0
                ground = compute_diagnostic_litter(ground, tile);
            else
                ground.CHILD = compute_diagnostic_CHILD(ground.CHILD, tile);
                ground = compute_diagnostic_litter(ground, tile);
            end
        end
        
        function ground = check_trigger(ground, tile)
            ground = check_trigger@GROUND_freezeC_RichardsEqW_Xice_seb_snow(ground, tile);
            
%             doy = tile.t - datenum(year(tile.t),1,1);
%             if doy > 152 && doy < 274 % Grow litter
%                 ground.STATVAR.layerThick(1) = ground.TEMP.target_layerThick + .05 .*(doy-152)./(274-152);
%                 ground.TEMP.litter = 1;
%             elseif doy >= 274
%                 ground.STATVAR.layerThick(1) = ground.TEMP.target_layerThick + .05 .*(1-(doy-274)./(366-274));
%                 ground.TEMP.litter = 1;
%             elseif doy < 1 && ground.TEMP.litter == 1
%                 ground.STATVAR.layerThick(1) = ground.TEMP.target_layerThick;
%             end            
%             ground = compute_diagnostic(ground, tile);
        end
        
        function z0 = get_z0_surface(ground)
            z0 = get_z0_surface@GROUND_freezeC_RichardsEqW_Xice_seb_snow(ground);
        end
        
        function albedo = get_albedo(ground)
            albedo = get_albedo@GROUND_freezeC_RichardsEqW_Xice_seb_snow(ground);
        end
        
        function Ts = get_surface_T(ground, tile)
            Ts = get_surface_T@GROUND_freezeC_RichardsEqW_Xice_seb_snow(ground, tile);
        end
        
        %----------
        
        function ground = compute_diagnostic_litter(ground, ~) %Identical to compute_diagnostic@GROUND_freezeC_RichardsEqW_Xice_seb, apart from thermCond
            %equilibrate water between matrix and Xwater within cells
            air = ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.XwaterIce - ground.STATVAR.waterIce - ground.STATVAR.mineral - ground.STATVAR.organic; 
            move_cells = (ground.STATVAR.Xwater > 0) & (air > 0);
            move_Xwater = min(ground.STATVAR.Xwater(move_cells), air(move_cells));
            ground.STATVAR.XwaterIce(move_cells) = max(0, ground.STATVAR.XwaterIce(move_cells) - move_Xwater);
            ground.STATVAR.waterIce(move_cells) = ground.STATVAR.waterIce(move_cells) + move_Xwater;
            ground.STATVAR.layerThick(move_cells) = ground.STATVAR.layerThick(move_cells) - move_Xwater ./  ground.STATVAR.area(move_cells);
            
            ground.STATVAR.layerThick = max(ground.STATVAR.layerThick, ...
                (ground.STATVAR.XwaterIce + ground.STATVAR.waterIce + ground.STATVAR.mineral + ground.STATVAR.organic) ./ ground.STATVAR.area);  %prevent rounding errors, would lead to wrong sign water fluxes in next prognostic step
            
            ground = get_T_water_freezeC_Xice(ground);
            
            ground = thermalConductivity_CLM4_5_Xice_litter(ground);
            ground = calculate_hydraulicConductivity_RichardsEq_Xice2(ground);
            
            ground = set_TEMP_2zero(ground);
        end
            
        %reset timestamp when changing TILES
        function ground = reset_timestamps(ground, tile)
            if ground.CHILD ~= 0
                ground.CHILD = reset_timestamps(ground.CHILD, tile);
            end
        end
        
    end
end