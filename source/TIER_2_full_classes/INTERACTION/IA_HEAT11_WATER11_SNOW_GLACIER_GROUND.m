%========================================================================
% CryoGrid INTERACTION (IA) class for heat conduction and water fluxes between a SNOW class
% and a GROUND class with bucket water
% contains function for SNOW CHILD phase 
% S. Westermann, October 2020
%========================================================================

classdef IA_HEAT11_WATER11_SNOW_GLACIER_GROUND < IA_WATER & IA_HEAT 
    
    methods
        
        function get_boundary_condition_m(ia_heat_water, tile)
            get_boundary_condition_HEAT_m(ia_heat_water);
            get_boundary_condition_BUCKET_SNOW_m(ia_heat_water);
            
            ia_heat_water.NEXT.TEMP.move_snow_yes_no = 1; 
        end
        
        %SNOW
        function get_IA_CHILD_boundary_condition_u(ia_heat_water, tile)
            get_boundary_condition_HEAT_IA_CHILD(ia_heat_water);
            get_boundary_condition_BUCKET_SNOW_m(ia_heat_water);
        end

        function move_snow2glacier(ia_heat_water, tile)
            move_cells = find((tile.t - ia_heat_water.PREVIOUS.STATVAR.time_snowfall) >  ia_heat_water.NEXT.PARA.min_snow_age .* 365.25);
            if ~isempty(move_cells)
                if move_cells(1,1) == 1 %always leave at least one cell in the SNOW class
                    move_cells(1,:) = [];
                end
                if ~isempty(move_cells)
                    snow = ia_heat_water.PREVIOUS;
                    glacier = ia_heat_water.NEXT;
                    
                    glacier.STATVAR.ice_fraction = [snow.STATVAR.ice(move_cells,1)./snow.STATVAR.layerThick(move_cells,1)./snow.STATVAR.area(move_cells,1); glacier.STATVAR.ice_fraction];
                    snow.STATVAR.target_density(move_cells,:) =  [];
                    snow.STATVAR.d(move_cells,:) = [];
                    snow.STATVAR.s(move_cells,:) = [];
                    snow.STATVAR.gs(move_cells,:) = [];
                    snow.STATVAR.time_snowfall(move_cells,:) = [];
                    snow.STATVAR.top_snow_date(move_cells,:) = [];
                    snow.STATVAR.bottom_snow_date(move_cells,:) = [];
                    
                    glacier.STATVAR.layerThick_glacier = [snow.STATVAR.layerThick(move_cells,1); glacier.STATVAR.layerThick_glacier];
                    glacier.STATVAR.layerThick_sediment = [zeros(size(move_cells,1),1); glacier.STATVAR.layerThick_sediment];
                    glacier.STATVAR.mineral = [zeros(size(move_cells,1),1); glacier.STATVAR.mineral];
                    glacier.STATVAR.organic = [zeros(size(move_cells,1),1); glacier.STATVAR.organic];
                    glacier.STATVAR.XwaterIce = [zeros(size(move_cells,1),1); glacier.STATVAR.XwaterIce];
                    glacier.STATVAR.Xwater = [zeros(size(move_cells,1),1); glacier.STATVAR.Xwater];
                    glacier.STATVAR.field_capacity_sediment = [zeros(size(move_cells,1),1) + glacier.STATVAR.field_capacity_sediment(1,1); glacier.STATVAR.field_capacity_sediment];
                    glacier.STATVAR.satHydraulicConductivity_sediment = [zeros(size(move_cells,1),1) + glacier.STATVAR.satHydraulicConductivity_sediment(1,1); glacier.STATVAR.satHydraulicConductivity_sediment];
                    glacier.TEMP.target_grid = [zeros(size(move_cells,1),1) + glacier.TEMP.target_grid(1,1); glacier.TEMP.target_grid];

                    common_variables = {'energy'; 'T'; 'water'; 'waterIce'; 'ice'; 'layerThick'; 'thermCond'; 'hydraulicConductivity'; 'area'};
                    for i=1:size(common_variables,1)
                        glacier.STATVAR.(common_variables{i,1}) = [snow.STATVAR.(common_variables{i,1})(move_cells,1); glacier.STATVAR.(common_variables{i,1})];
                        snow.STATVAR.(common_variables{i,1})(move_cells,:) =  [];
                    end
                    snow = compute_diagnostic(snow, tile);
                    glacier = compute_diagnostic(glacier, tile);
                end
            end
        end
        
        function remove_excessWater_CHILD(ia_heat_water) %move excessWater from SNOW to water and Xwater from PARENT, 0 energy transfer since meltwater must be zero
            space_left = max(0,ia_heat_water.NEXT.STATVAR.layerThick(1) .* ia_heat_water.NEXT.STATVAR.area(1) - ia_heat_water.NEXT.STATVAR.mineral(1) ...
                - ia_heat_water.NEXT.STATVAR.organic(1) - ia_heat_water.NEXT.STATVAR.waterIce(1)); 
            water_in = min(space_left, ia_heat_water.PREVIOUS.STATVAR.excessWater);
            Xwater_in = max(0, ia_heat_water.PREVIOUS.STATVAR.excessWater - water_in);
            
            ia_heat_water.NEXT.STATVAR.waterIce(1) = ia_heat_water.NEXT.STATVAR.waterIce(1) + water_in + Xwater_in;
            ia_heat_water.NEXT.STATVAR.XwaterIce(1) = ia_heat_water.NEXT.STATVAR.XwaterIce(1) + Xwater_in;
            ia_heat_water.NEXT.STATVAR.Xwater(1) = ia_heat_water.NEXT.STATVAR.Xwater(1) + Xwater_in;
            ia_heat_water.NEXT.STATVAR.layerThick(1) = ia_heat_water.NEXT.STATVAR.layerThick(1) + Xwater_in ./ ia_heat_water.NEXT.STATVAR.area(1);
            ia_heat_water.NEXT.STATVAR.energy(1,1) = ia_heat_water.NEXT.STATVAR.energy(1,1) + ia_heat_water.PREVIOUS.STATVAR.excessWater_energy;
            ia_heat_water.PREVIOUS.STATVAR.excessWater = 0;
            ia_heat_water.PREVIOUS.STATVAR.excessWater_energy = 0;
        end
    end
end